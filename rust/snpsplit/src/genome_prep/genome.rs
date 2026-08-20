//! Reading the reference genome, applying SNPs, and writing the modified chromosomes.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use anyhow::{Context, Result};

/// One SNP as the per-chromosome track records it: position, reference base, SNP base.
pub struct Snp {
    pub pos: usize,
    pub reference: u8,
    pub alternative: u8,
}

/// The reference genome, chromosome name to uppercased sequence.
pub type Genome = BTreeMap<String, Vec<u8>>;

/// The chromosome name a FastA header declares.
///
/// Bowtie takes the first whitespace-delimited token after `>`, and so does this, so that a
/// header carrying a description does not produce a chromosome named after the description.
pub fn extract_chromosome_name(header: &str) -> Result<String> {
    let Some(rest) = header.strip_prefix('>') else {
        anyhow::bail!(
            "The specified chromosome ({header}) file doesn't seem to be in FASTA format as required!\n"
        );
    };
    Ok(rest.split_whitespace().next().unwrap_or("").to_string())
}

/// Read every FastA file in `folder` into memory.
///
/// `.fa` first, falling back to `.fasta`, matching the Perl's two globs. Files are taken in
/// sorted order, which is what the shell glob the Perl uses produces.
pub fn read_genome_into_memory(
    folder: &Path,
    out: &mut impl Write,
    err: &mut impl Write,
) -> Result<Genome> {
    writeln!(
        err,
        "Now reading in and storing sequence information of the genome specified in: {}\n",
        folder.display()
    )?;

    let mut names = fasta_files(folder, "fa")?;
    if names.is_empty() {
        names = fasta_files(folder, "fasta")?;
    }
    if names.is_empty() {
        anyhow::bail!(
            "The specified reference genome folder {} does not contain any sequence files in FastA format (with .fa or .fasta file extensions)\n",
            folder.display()
        );
    }

    let mut genome: Genome = BTreeMap::new();

    for name in names {
        let path = folder.join(&name);
        let file = File::open(&path)
            .with_context(|| format!("Failed to read from sequence file {name} "))?;
        let mut lines = BufReader::new(file).lines();

        let first = lines.next().transpose()?.unwrap_or_default();
        let first = first.trim_end_matches('\r');
        let mut chromosome_name = extract_chromosome_name(first)?;
        let mut sequence: Vec<u8> = Vec::new();

        for line in lines {
            let line = line?;
            let line = line.trim_end_matches('\r');
            if line.starts_with('>') {
                store(&mut genome, &chromosome_name, &sequence, &name, out, true)?;
                sequence = Vec::new();
                chromosome_name = extract_chromosome_name(line)?;
            } else {
                sequence.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
            }
        }

        store(&mut genome, &chromosome_name, &sequence, &name, out, false)?;
    }

    writeln!(out)?;
    Ok(genome)
}

/// Store one chromosome, refusing a duplicate name.
///
/// `mid_file` distinguishes the two call sites in the Perl, which print the same line with a
/// different terminator before dying. The difference is observable, so it is kept.
fn store(
    genome: &mut Genome,
    name: &str,
    sequence: &[u8],
    filename: &str,
    out: &mut impl Write,
    mid_file: bool,
) -> Result<()> {
    if genome.contains_key(name) {
        if mid_file {
            writeln!(out, "chr {name} ({} bp)", sequence.len())?;
        } else {
            write!(out, "chr {name} ({} bp)\t", sequence.len())?;
        }
        anyhow::bail!(
            "Exiting because chromosome name already exists. Please make sure all chromosomes have a unique name{}\n",
            if mid_file { "!" } else { "." }
        );
    }

    if sequence.is_empty() {
        let where_ = if mid_file {
            format!("multi-fasta file {filename}")
        } else {
            format!("file {filename}")
        };
        // Written to stderr by the Perl, and this function only has stdout, so the caller's
        // ordering is preserved by emitting it directly.
        eprintln!("Chromosome {name} in the {where_} did not contain any sequence information!");
    }

    writeln!(out, "chr {name} ({} bp)", sequence.len())?;
    genome.insert(name.to_string(), sequence.to_vec());
    Ok(())
}

fn fasta_files(folder: &Path, extension: &str) -> Result<Vec<String>> {
    let mut names: Vec<String> = std::fs::read_dir(folder)
        .with_context(|| format!("Failed to read genome folder {}", folder.display()))?
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.file_name().to_string_lossy().to_string())
        .filter(|name| name.ends_with(&format!(".{extension}")))
        .collect();
    names.sort();
    Ok(names)
}

/// Read the per-chromosome SNP track written by the filtering step.
///
/// Returns an empty list when the file is absent: some chromosomes legitimately carry no
/// SNPs, and the sequence is then written back out unmodified.
pub fn read_snps(parent: &Path, chr: &str, strain: &str, err: &mut impl Write) -> Result<Vec<Snp>> {
    let folder = parent.join(format!("SNPs_{strain}"));
    if !folder.is_dir() {
        anyhow::bail!(
            "Folder >>{}<< doesn't exist. Try losing the option --skip_filtering to generate the folder and SNP files from the VCF file\n\n",
            folder.display()
        );
    }

    let file = folder.join(format!("chr{chr}.txt"));
    if !file.exists() {
        writeln!(
            err,
            "Couldn't find SNP file for chromosome '{chr}' '{}' didn't exist. Skipping...",
            file.display()
        )?;
        return Ok(Vec::new());
    }
    writeln!(err, "Reading SNPs from file {}", file.display())?;

    let mut snps = Vec::new();
    for line in BufReader::new(File::open(&file)?).lines() {
        let line = line?;
        let line = line.replace('\r', "");
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let (Some(pos), Some(strand), Some(allele)) = (fields.get(2), fields.get(3), fields.get(4))
        else {
            continue;
        };
        if allele.is_empty() {
            continue;
        }

        let Some((reference, alternative)) = parse_allele(allele) else {
            writeln!(
                err,
                "Skipping allele '{allele}' as it appears to contain non DNA bases (only G,A,T,C allowed)"
            )?;
            continue;
        };

        let (reference, alternative) = if *strand == "-1" {
            (complement(reference), complement(alternative))
        } else {
            (reference, alternative)
        };

        let Ok(pos) = pos.parse::<usize>() else {
            continue;
        };
        snps.push(Snp {
            pos,
            reference,
            alternative,
        });
    }

    snps.sort_by_key(|s| s.pos);
    Ok(snps)
}

/// Split a `X/Y` allele, rejecting anything that is not a pair of unambiguous bases.
pub fn parse_allele(allele: &str) -> Option<(u8, u8)> {
    let bytes = allele.as_bytes();
    if bytes.len() != 3 || bytes[1] != b'/' {
        return None;
    }
    let ok = |b: u8| matches!(b, b'G' | b'A' | b'T' | b'C');
    (ok(bytes[0]) && ok(bytes[2])).then(|| (bytes[0], bytes[2]))
}

/// Complement a single base, leaving anything else alone, as Perl's `tr/GATC/CTAG/` does.
pub fn complement(base: u8) -> u8 {
    match base {
        b'G' => b'C',
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        other => other,
    }
}

/// Counts produced by applying one chromosome's SNPs.
#[derive(Default)]
pub struct Applied {
    pub total: usize,
    pub already_carried: usize,
    pub mismatched: usize,
    pub new_n: usize,
    pub new_snp: usize,
}

/// Apply `snps` to `sequence`, producing the N-masked and full-sequence versions asked for.
///
/// Positions are 1-based. A SNP whose reference base disagrees with the genome is skipped and
/// counted: it means the VCF and the FastA are not the same assembly, which is worth knowing.
pub fn apply_snps(
    sequence: &[u8],
    snps: &[Snp],
    nmasking: bool,
    full_sequence: bool,
) -> (Applied, Option<Vec<u8>>, Option<Vec<u8>>) {
    let mut counts = Applied::default();
    let mut full = sequence.to_vec();
    let mut masked = if nmasking {
        Some(sequence.to_vec())
    } else {
        None
    };

    let mut last_pos = 0usize;
    for snp in snps {
        counts.total += 1;
        if snp.pos == last_pos {
            continue;
        }
        last_pos = snp.pos;

        let Some(&base) = full.get(snp.pos.wrapping_sub(1)) else {
            counts.mismatched += 1;
            continue;
        };

        if base == snp.alternative {
            counts.already_carried += 1;
            continue;
        }
        if base != snp.reference {
            counts.mismatched += 1;
            continue;
        }

        if let Some(masked) = masked.as_mut() {
            masked[snp.pos - 1] = b'N';
            counts.new_n += 1;
        }
        if full_sequence {
            full[snp.pos - 1] = snp.alternative;
            counts.new_snp += 1;
        }
    }

    let full = full_sequence.then_some(full);
    (counts, masked, full)
}

/// Write one chromosome out at 100 bases per line.
pub fn write_chromosome(
    parent: &Path,
    chr: &str,
    sequence: &[u8],
    nmasked: bool,
    strain: &str,
    err: &mut impl Write,
) -> Result<()> {
    if nmasked {
        writeln!(err, "Writing modified chromosome (N-masking)")?;
    } else {
        writeln!(err, "Writing modified chromosome (incorporating SNPs)")?;
    }

    let (kind, outfile) = if nmasked {
        ("N-masked", format!("chr{chr}.N-masked.fa"))
    } else {
        ("full_sequence", format!("chr{chr}.SNPs_introduced.fa"))
    };

    let folder = parent.join(format!("{strain}_{kind}"));
    let path = folder.join(&outfile);
    if nmasked {
        writeln!(err, "Writing N-masked output to: {}", path.display())?;
    } else {
        writeln!(err, "Writing full sequence output to: {}", path.display())?;
    }

    if !folder.is_dir() {
        let _ = std::fs::create_dir(&folder);
    }

    let mut fh = std::io::BufWriter::new(
        File::create(&path)
            .with_context(|| format!("Failed to write to file {}: ", path.display()))?,
    );
    writeln!(fh, ">{chr}")?;

    // The Perl's loop condition leaves the final partial line to a separate write, and a
    // sequence whose length is an exact multiple of 100 therefore ends with a full line
    // rather than an empty one. An empty sequence writes a single blank line.
    let mut pos = 0usize;
    while pos + 100 < sequence.len() {
        fh.write_all(&sequence[pos..pos + 100])?;
        fh.write_all(b"\n")?;
        pos += 100;
    }
    fh.write_all(&sequence[pos..])?;
    fh.write_all(b"\n")?;
    fh.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_header_name_stops_at_whitespace() {
        assert_eq!(extract_chromosome_name(">1").unwrap(), "1");
        assert_eq!(
            extract_chromosome_name(">chr1 some description").unwrap(),
            "chr1"
        );
    }

    #[test]
    fn a_line_that_is_not_a_header_is_rejected() {
        assert!(extract_chromosome_name("ACGT").is_err());
    }

    #[test]
    fn only_a_pair_of_unambiguous_bases_is_an_allele() {
        assert_eq!(parse_allele("T/A"), Some((b'T', b'A')));
        assert_eq!(parse_allele("T/N"), None);
        assert_eq!(parse_allele("TT/A"), None);
        assert_eq!(parse_allele("T|A"), None);
        assert_eq!(parse_allele(""), None);
    }

    #[test]
    fn complementing_leaves_anything_that_is_not_a_base_alone() {
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn a_snp_whose_reference_base_disagrees_is_skipped_and_counted() {
        let seq = b"ACGT".to_vec();
        let snps = vec![Snp {
            pos: 1,
            reference: b'G',
            alternative: b'T',
        }];
        let (counts, masked, _) = apply_snps(&seq, &snps, true, false);
        assert_eq!(counts.mismatched, 1);
        assert_eq!(counts.new_n, 0);
        assert_eq!(masked.unwrap(), seq);
    }

    #[test]
    fn a_position_already_carrying_the_snp_base_is_left_alone() {
        let seq = b"ACGT".to_vec();
        let snps = vec![Snp {
            pos: 1,
            reference: b'G',
            alternative: b'A',
        }];
        let (counts, _, _) = apply_snps(&seq, &snps, true, false);
        assert_eq!(counts.already_carried, 1);
        assert_eq!(counts.new_n, 0);
    }

    /// Duplicate positions are skipped but still counted in the total, which is why the total
    /// does not always reconcile with the applied and skipped counts. Reproduced, and raised
    /// upstream separately.
    #[test]
    fn a_duplicate_position_is_counted_but_applied_once() {
        let seq = b"ACGT".to_vec();
        let snps = vec![
            Snp {
                pos: 1,
                reference: b'A',
                alternative: b'C',
            },
            Snp {
                pos: 1,
                reference: b'A',
                alternative: b'C',
            },
        ];
        let (counts, _, _) = apply_snps(&seq, &snps, true, false);
        assert_eq!(counts.total, 2);
        assert_eq!(counts.new_n, 1);
    }
}
