//! Dual hybrid genomes: comparing two strains and building a genome whose reference is
//! strain 1 and whose SNPs are strain 2.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};

use super::genome::{Snp, complement, parse_allele};
use super::vcf::Blacklisted;

/// One strain-1 SNP, and whether strain 2 also carried it.
struct Position {
    reference: u8,
    snp: u8,
    seen_in_strain2: usize,
}

/// Counts the comparison reports.
#[derive(Default)]
pub struct Comparison {
    pub strain_positions: usize,
    pub strain2_positions: usize,
    pub same: usize,
    pub different: usize,
    pub confidence_discrepancy: usize,
    pub unique_ref: usize,
    pub unique_snp: usize,
}

fn open_maybe_gzipped(path: &Path) -> Result<Box<dyn BufRead>> {
    let file =
        File::open(path).with_context(|| format!("Failed to read from file {}", path.display()))?;
    let reader: Box<dyn Read> = if path.to_string_lossy().ends_with(".gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

/// Compare the two strains' SNP lists and write the four dual hybrid annotation files.
///
/// Returns the name of the combined annotation, which the third genome is built from.
#[allow(clippy::too_many_arguments)]
pub fn determine_snps_between_strains(
    strain: &str,
    strain2: &str,
    genome_build: &str,
    snp_file_strain: &Path,
    snp_file_strain2: &Path,
    blacklist: &BTreeMap<String, Blacklisted>,
    report: &mut impl Write,
    err: &mut impl Write,
) -> Result<(String, Comparison)> {
    writeln!(
        err,
        "Determining new Ref [{strain}] and SNP [{strain2}] annotations"
    )?;
    writeln!(err, "{}", "=".repeat(60))?;
    writeln!(err)?;

    let out_strain_name = format!("{strain}_specific_SNPs.{genome_build}.txt");
    writeln!(
        err,
        "Writing {strain} specific SNPs (relative to the {genome_build} reference) to >>{out_strain_name}<<"
    )?;
    let out_strain2_name = format!("{strain2}_specific_SNPs.{genome_build}.txt");
    writeln!(
        err,
        "Writing {strain2} specific SNPs (relative to the {genome_build} reference) to >>{out_strain2_name}<<"
    )?;
    let out_common_name = format!("{strain}_{strain2}_SNPs_in_common.{genome_build}.txt");
    writeln!(
        err,
        "Writing SNPs in common between {strain} and {strain2} (relative to the {genome_build} reference) to >>{out_common_name}<<"
    )?;
    let all_name = format!("all_{strain2}_SNPs_{strain}_reference.based_on_{genome_build}.txt");
    writeln!(
        err,
        "Writing all new SNPs >>{strain}/{strain2} to >>{all_name}<<\n"
    )?;

    if !snp_file_strain.exists() {
        anyhow::bail!(
            "Expected SNP file [{}] for strain {strain} did not exist! Please make sure that it is present in the current working directory\n\n",
            snp_file_strain.display()
        );
    }
    if !snp_file_strain2.exists() {
        // The stray brace is the Perl's, in an interpolated file name. Kept because the
        // message is diffed.
        anyhow::bail!(
            "Expected SNP file 2 [{}}}] for strain {strain2} did not exist! Please make sure that it is present in the current working directory\n\n",
            snp_file_strain2.display()
        );
    }

    let mut out_strain = File::create(&out_strain_name)?;
    let mut out_strain2 = File::create(&out_strain2_name)?;
    let mut out_common = File::create(&out_common_name)?;
    let mut all = File::create(&all_name)?;

    let mut counts = Comparison::default();
    // chr -> pos -> position. BTreeMap so the final sweep is in the sorted order the Perl
    // takes with `sort keys` and `sort {$a <=> $b}`.
    let mut snps: BTreeMap<String, BTreeMap<u64, Position>> = BTreeMap::new();

    writeln!(
        err,
        "Storing SNP positions for strain {strain} provided in '{}'",
        snp_file_strain.display()
    )?;
    for line in open_maybe_gzipped(snp_file_strain)?.lines() {
        let line = line?;
        counts.strain_positions += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(chr), Some(pos), Some(diff)) = (fields.get(1), fields.get(2), fields.get(4))
        else {
            continue;
        };
        let Some((reference, snp)) = split_pair(diff) else {
            continue;
        };
        let Ok(pos) = pos.parse::<u64>() else {
            continue;
        };
        snps.entry(chr.to_string()).or_default().insert(
            pos,
            Position {
                reference,
                snp,
                seen_in_strain2: 0,
            },
        );
    }
    writeln!(
        err,
        "Stored {} positions in total\n",
        counts.strain_positions
    )?;

    writeln!(
        err,
        "Now reading and comparing SNP positions for strain {strain2} provided in '{}'",
        snp_file_strain2.display()
    )?;
    for line in open_maybe_gzipped(snp_file_strain2)?.lines() {
        let line = line?;
        counts.strain2_positions += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(chr), Some(pos_text), Some(diff)) = (fields.get(1), fields.get(2), fields.get(4))
        else {
            continue;
        };
        let Some((reference, snp)) = split_pair(diff) else {
            continue;
        };
        let Ok(pos) = pos_text.parse::<u64>() else {
            continue;
        };
        let location = format!("{chr}:{pos_text}");

        match snps.get_mut(*chr).and_then(|m| m.get_mut(&pos)) {
            Some(existing) => {
                existing.seen_in_strain2 += 1;
                if reference != existing.reference {
                    writeln!(err, "reference was different for the same position!!!")?;
                }
                if snp == existing.snp {
                    counts.same += 1;
                    writeln!(out_common, "{line}")?;
                } else {
                    counts.different += 1;
                    writeln!(
                        all,
                        "{}\t{chr}\t{pos_text}\t1\t{}/{}",
                        counts.different, existing.snp as char, snp as char
                    )?;
                }
            }
            None => {
                // Absent from strain 1's high-confidence list is not the same as absent from
                // the VCF: a position that was homozygous but low confidence in strain 1 is
                // dropped rather than treated as strain-2 specific.
                if let Some(entry) = blacklist.get(&location)
                    && let Some(filter) = &entry.strain1_filter
                    && filter != "1"
                {
                    counts.confidence_discrepancy += 1;
                    continue;
                }
                counts.unique_snp += 1;
                writeln!(out_strain2, "{line}")?;
                writeln!(all, "{line}")?;
            }
        }
    }

    writeln!(
        err,
        "Finally, looking at new reference [{strain}] specific reads..."
    )?;

    for (chr, positions) in &snps {
        for (pos, position) in positions {
            let location = format!("{chr}:{pos}");
            if position.seen_in_strain2 == 0 {
                if let Some(entry) = blacklist.get(&location)
                    && let Some(filter) = &entry.strain2_filter
                    && filter != "1"
                {
                    counts.confidence_discrepancy += 1;
                    continue;
                }
                counts.unique_ref += 1;
                let line = format!(
                    "{strain}_{}\t{chr}\t{pos}\t1\t{}/{}",
                    counts.unique_ref, position.snp as char, position.reference as char
                );
                writeln!(out_strain, "{line}")?;
                writeln!(all, "{line}")?;
            }
            if position.seen_in_strain2 >= 2 {
                anyhow::bail!(
                    "SNP was present at least twice: {chr}\t{pos}\tcount: {}\n\n",
                    position.snp as char
                );
            }
        }
    }

    let summary = format!(
        "Looked at positions from new Reference strain [{strain}]:\t\t{}\nCompared positions from new SNP strain [{strain2}]:\t\t{}\n{}\nSNPs were the same in Ref and SNP genome (not written out):\t{}\nSNPs were present in both Ref and SNP genome but had a different sequence:\t{}\nSNPs were low confidence in one strain and thus ignored:\t{}\nSNPs were unique to Ref [{strain}]:\t\t\t\t{}\nSNPs were unique to SNP [{strain2}]:\t\t\t\t{}\n\n",
        counts.strain_positions,
        counts.strain2_positions,
        "=".repeat(54),
        counts.same,
        counts.different,
        counts.confidence_discrepancy,
        counts.unique_ref,
        counts.unique_snp,
    );
    write!(err, "\n{summary}")?;
    write!(report, "{summary}")?;

    Ok((all_name, counts))
}

fn split_pair(diff: &str) -> Option<(u8, u8)> {
    let (a, b) = diff.split_once('/')?;
    Some((*a.as_bytes().first()?, *b.as_bytes().first()?))
}

/// Read the combined strain1/strain2 annotation back in, for the third genome.
pub fn read_new_snp_annotation(
    path: &Path,
    strain: &str,
    strain2: &str,
    err: &mut impl Write,
) -> Result<BTreeMap<String, Vec<Snp>>> {
    writeln!(
        err,
        "Reading {strain}/{strain2} SNPs from file '{}'",
        path.display()
    )?;
    if !path.exists() {
        anyhow::bail!("Couldn't find SNP file '{}'\n", path.display());
    }

    let mut by_chr: BTreeMap<String, BTreeMap<usize, Snp>> = BTreeMap::new();
    for line in open_maybe_gzipped(path)?.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(chr), Some(pos), Some(strand), Some(allele)) =
            (fields.get(1), fields.get(2), fields.get(3), fields.get(4))
        else {
            continue;
        };
        let Some((reference, snp)) = parse_allele(allele) else {
            writeln!(err, "Skipping allele '{allele}'")?;
            continue;
        };
        let (reference, snp) = if *strand == "-1" {
            (complement(reference), complement(snp))
        } else {
            (reference, snp)
        };
        let Ok(pos) = pos.parse::<usize>() else {
            continue;
        };
        by_chr.entry(chr.to_string()).or_default().insert(
            pos,
            Snp {
                pos,
                reference,
                alternative: snp,
            },
        );
    }

    Ok(by_chr
        .into_iter()
        .map(|(chr, positions)| (chr, positions.into_values().collect()))
        .collect())
}
