//! The three sorting modes and the reports they produce.

use std::io::Write;

use anyhow::Result;
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

use crate::io::RecordReader;

use super::outputs::Outputs;
use super::{Allele, allele_of, percentage};

fn name_of(record: &RecordBuf) -> Vec<u8> {
    record
        .name()
        .map(|n| {
            let bytes: &[u8] = n.as_ref();
            bytes.to_vec()
        })
        .unwrap_or_default()
}

fn require_allele(record: &RecordBuf) -> Result<Allele> {
    allele_of(record).ok_or_else(|| {
        anyhow::anyhow!(
            "Failed to extract XX:Z tag from line:\n{:?}\n",
            name_of(record)
        )
    })
}

/// What a run counted, in the shape its mode reports.
pub enum Counts {
    SingleEnd {
        total: usize,
        unassigned: usize,
        genome1: usize,
        genome2: usize,
        conflicting: usize,
    },
    PairedEnd {
        total: usize,
        pairs: usize,
        singletons: usize,
        unassigned: usize,
        unassigned_pairs: usize,
        unassigned_singletons: usize,
        genome1: usize,
        g1_pairs: usize,
        g1_singletons: usize,
        genome2: usize,
        g2_pairs: usize,
        g2_singletons: usize,
        conflicting: usize,
        conflicting_pairs: usize,
        conflicting_singletons: usize,
    },
    HiC {
        total: usize,
        unassigned: usize,
        genome1: usize,
        genome2: usize,
        conflicting: usize,
        g1_ua: usize,
        ua_g1: usize,
        g2_ua: usize,
        ua_g2: usize,
        g1_g2: usize,
        g2_g1: usize,
    },
}

/// Single-end: every alignment is classified on its own.
pub fn process_single_end(
    reader: RecordReader,
    header: &Header,
    outputs: &mut Outputs,
    verbose: bool,
    err: &mut impl Write,
) -> Result<Counts> {
    let (mut unassigned, mut genome1, mut genome2, mut conflicting) = (0, 0, 0, 0);
    let mut total = 0usize;

    for record in reader {
        let record = record?;
        total += 1;
        if total.is_multiple_of(1_000_000) {
            writeln!(err, "Processed {total} lines so far")?;
        }
        if verbose {
            println!("{}", crate::io::render_sam_line(header, &record)?);
        }

        match require_allele(&record)? {
            Allele::Unassigned => {
                unassigned += 1;
                Outputs::write(&mut outputs.unassigned, &record)?;
            }
            Allele::Genome1 => {
                genome1 += 1;
                Outputs::write(&mut outputs.genome1, &record)?;
            }
            Allele::Genome2 => {
                genome2 += 1;
                Outputs::write(&mut outputs.genome2, &record)?;
            }
            Allele::Conflicting => {
                conflicting += 1;
                if let Some(sink) = outputs.conflicting.as_mut() {
                    Outputs::write(sink, &record)?;
                }
            }
        }
    }

    Ok(Counts::SingleEnd {
        total,
        unassigned,
        genome1,
        genome2,
        conflicting,
    })
}

/// Paired-end: consecutive alignments sharing a read name are a pair, anything else is a
/// singleton. The input is name-sorted, so adjacency is the whole test.
pub fn process_paired_end(
    reader: RecordReader,
    _header: &Header,
    outputs: &mut Outputs,
    singletons_separately: bool,
    err: &mut impl Write,
) -> Result<Counts> {
    let mut c = PairedCounts::default();
    let mut held: Option<RecordBuf> = None;

    for record in reader {
        let record = record?;
        match held.take() {
            None => held = Some(record),
            Some(previous) => {
                if name_of(&previous) == name_of(&record) {
                    c.total += 1;
                    c.pairs += 1;
                    let a = require_allele(&previous)?;
                    let b = require_allele(&record)?;
                    classify_pair(a, b, &previous, &record, outputs, &mut c)?;
                    held = None;
                } else {
                    c.total += 1;
                    c.singletons += 1;
                    let allele = require_allele(&previous)?;
                    classify_singleton(allele, &previous, outputs, singletons_separately, &mut c)?;
                    held = Some(record);
                }
                if c.total.is_multiple_of(1_000_000) {
                    writeln!(err, "Processed {} lines so far", c.total)?;
                }
            }
        }
    }

    match held {
        Some(last) => {
            c.total += 1;
            c.singletons += 1;
            let allele = require_allele(&last)?;
            classify_singleton(allele, &last, outputs, singletons_separately, &mut c)?;
        }
        None => writeln!(
            err,
            "Last read was a read pair which has already been processed. all done\n"
        )?,
    }

    Ok(Counts::PairedEnd {
        total: c.total,
        pairs: c.pairs,
        singletons: c.singletons,
        unassigned: c.unassigned,
        unassigned_pairs: c.unassigned_pairs,
        unassigned_singletons: c.unassigned_singletons,
        genome1: c.genome1,
        g1_pairs: c.g1_pairs,
        g1_singletons: c.g1_singletons,
        genome2: c.genome2,
        g2_pairs: c.g2_pairs,
        g2_singletons: c.g2_singletons,
        conflicting: c.conflicting,
        conflicting_pairs: c.conflicting_pairs,
        conflicting_singletons: c.conflicting_singletons,
    })
}

#[derive(Default)]
struct PairedCounts {
    total: usize,
    pairs: usize,
    singletons: usize,
    unassigned: usize,
    unassigned_pairs: usize,
    unassigned_singletons: usize,
    genome1: usize,
    g1_pairs: usize,
    g1_singletons: usize,
    genome2: usize,
    g2_pairs: usize,
    g2_singletons: usize,
    conflicting: usize,
    conflicting_pairs: usize,
    conflicting_singletons: usize,
}

/// A pair goes wherever its more informative mate points, and a G1/G2 disagreement is
/// conflicting rather than a mix.
fn classify_pair(
    a: Allele,
    b: Allele,
    first: &RecordBuf,
    second: &RecordBuf,
    outputs: &mut Outputs,
    c: &mut PairedCounts,
) -> Result<()> {
    use Allele::*;
    match (a, b) {
        (Unassigned, Unassigned) => {
            c.unassigned += 1;
            c.unassigned_pairs += 1;
            Outputs::write(&mut outputs.unassigned, first)?;
            Outputs::write(&mut outputs.unassigned, second)?;
        }
        (Genome1, Unassigned) | (Unassigned, Genome1) | (Genome1, Genome1) => {
            c.genome1 += 1;
            c.g1_pairs += 1;
            Outputs::write(&mut outputs.genome1, first)?;
            Outputs::write(&mut outputs.genome1, second)?;
        }
        (Genome2, Unassigned) | (Unassigned, Genome2) | (Genome2, Genome2) => {
            c.genome2 += 1;
            c.g2_pairs += 1;
            Outputs::write(&mut outputs.genome2, first)?;
            Outputs::write(&mut outputs.genome2, second)?;
        }
        _ => {
            // Both the G1/G2 disagreement and any CF mate land here.
            c.conflicting += 1;
            c.conflicting_pairs += 1;
            if let Some(sink) = outputs.conflicting.as_mut() {
                Outputs::write(sink, first)?;
                Outputs::write(sink, second)?;
            }
        }
    }
    Ok(())
}

fn classify_singleton(
    allele: Allele,
    record: &RecordBuf,
    outputs: &mut Outputs,
    separately: bool,
    c: &mut PairedCounts,
) -> Result<()> {
    match allele {
        Allele::Unassigned => {
            c.unassigned += 1;
            c.unassigned_singletons += 1;
            let sink = if separately {
                outputs
                    .unassigned_st
                    .as_mut()
                    .expect("opened with --singletons")
            } else {
                &mut outputs.unassigned
            };
            Outputs::write(sink, record)?;
        }
        Allele::Genome1 => {
            c.genome1 += 1;
            c.g1_singletons += 1;
            let sink = if separately {
                outputs
                    .genome1_st
                    .as_mut()
                    .expect("opened with --singletons")
            } else {
                &mut outputs.genome1
            };
            Outputs::write(sink, record)?;
        }
        Allele::Genome2 => {
            c.genome2 += 1;
            c.g2_singletons += 1;
            let sink = if separately {
                outputs
                    .genome2_st
                    .as_mut()
                    .expect("opened with --singletons")
            } else {
                &mut outputs.genome2
            };
            Outputs::write(sink, record)?;
        }
        Allele::Conflicting => {
            c.conflicting += 1;
            c.conflicting_singletons += 1;
            let sink = if separately {
                outputs.conflicting_st.as_mut()
            } else {
                outputs.conflicting.as_mut()
            };
            if let Some(sink) = sink {
                Outputs::write(sink, record)?;
            }
        }
    }
    Ok(())
}

/// Hi-C: paired by definition, and every combination of the two tags is its own destination.
pub fn process_hic(
    reader: RecordReader,
    _header: &Header,
    outputs: &mut Outputs,
    err: &mut impl Write,
) -> Result<Counts> {
    let (mut unassigned, mut genome1, mut genome2, mut conflicting) = (0, 0, 0, 0);
    let (mut g1_ua, mut ua_g1, mut g2_ua, mut ua_g2, mut g1_g2, mut g2_g1) = (0, 0, 0, 0, 0, 0);
    let mut total = 0usize;
    let mut held: Option<RecordBuf> = None;

    for record in reader {
        let record = record?;
        let Some(previous) = held.take() else {
            held = Some(record);
            continue;
        };

        if name_of(&previous) != name_of(&record) {
            anyhow::bail!(
                "Hi-C data has to be paired-end by definition, however the reads did not have the same read ID:\nRead 1: {:?}\nRead 2: {:?}\n\n",
                name_of(&previous),
                name_of(&record)
            );
        }

        total += 1;
        let a = require_allele(&previous)?;
        let b = require_allele(&record)?;

        use Allele::*;
        let sink = match (a, b) {
            (Unassigned, Unassigned) => {
                unassigned += 1;
                Some(&mut outputs.unassigned)
            }
            (Genome1, Genome1) => {
                genome1 += 1;
                Some(&mut outputs.genome1)
            }
            (Genome2, Genome2) => {
                genome2 += 1;
                Some(&mut outputs.genome2)
            }
            // Checked before the mixed combinations, as in the Perl: a CF on either mate
            // wins over what the other mate says.
            (Conflicting, _) | (_, Conflicting) => {
                conflicting += 1;
                outputs.conflicting.as_mut()
            }
            (Genome1, Unassigned) => {
                g1_ua += 1;
                outputs.g1_ua.as_mut()
            }
            (Unassigned, Genome1) => {
                ua_g1 += 1;
                outputs.g1_ua.as_mut()
            }
            (Genome2, Unassigned) => {
                g2_ua += 1;
                outputs.g2_ua.as_mut()
            }
            (Unassigned, Genome2) => {
                ua_g2 += 1;
                outputs.g2_ua.as_mut()
            }
            (Genome1, Genome2) => {
                g1_g2 += 1;
                outputs.g1_g2.as_mut()
            }
            (Genome2, Genome1) => {
                g2_g1 += 1;
                outputs.g1_g2.as_mut()
            }
        };

        if let Some(sink) = sink {
            Outputs::write(sink, &previous)?;
            Outputs::write(sink, &record)?;
        }

        if total.is_multiple_of(1_000_000) {
            writeln!(err, "Processed {total} lines so far")?;
        }
    }

    if held.is_some() {
        writeln!(err, "Last read was a singleton, skipping...")?;
    }

    Ok(Counts::HiC {
        total,
        unassigned,
        genome1,
        genome2,
        conflicting,
        g1_ua,
        ua_g1,
        g2_ua,
        ua_g2,
        g1_g2,
        g2_g1,
    })
}

impl Counts {
    /// Write the summary to stderr, the report file, and the YAML file.
    pub fn write(&self, outputs: &mut Outputs, err: &mut impl Write) -> Result<()> {
        match self {
            Counts::SingleEnd {
                total,
                unassigned,
                genome1,
                genome2,
                conflicting,
            } => {
                let p = |n: &usize| percentage(*n, *total);
                let body = format!(
                    "\n\nAllele-specific single-end sorting report\n{}\nRead alignments processed in total:\t\t{total}\nReads were unassignable:\t\t\t{unassigned} ({}%)\nReads were specific for genome 1:\t\t{genome1} ({}%)\nReads were specific for genome 2:\t\t{genome2} ({}%)\nReads contained conflicting SNP information:\t{conflicting} ({})\n\n\n",
                    "=".repeat(41),
                    p(unassigned),
                    p(genome1),
                    p(genome2),
                    // No percent sign here, unlike every other line. The Perl omits it and
                    // the fixtures diff it.
                    p(conflicting),
                );
                write!(err, "{body}")?;
                write!(outputs.report, "{body}")?;

                outputs.yaml("SE_total_reads", total);
                outputs.yaml("SE_unassignable", unassigned);
                outputs.yaml("SE_percent_unassignable", p(unassigned));
                outputs.yaml("SE_genome1", genome1);
                outputs.yaml("SE_percent_genome1", p(genome1));
                outputs.yaml("SE_genome2", genome2);
                outputs.yaml("SE_percent_genome2", p(genome2));
                outputs.yaml("SE_conflicting", conflicting);
                outputs.yaml("SE_percent_conflicting", p(conflicting));
            }

            Counts::PairedEnd {
                total,
                pairs,
                singletons,
                unassigned,
                unassigned_pairs,
                unassigned_singletons,
                genome1,
                g1_pairs,
                g1_singletons,
                genome2,
                g2_pairs,
                g2_singletons,
                conflicting,
                conflicting_pairs,
                conflicting_singletons,
            } => {
                let p = |n: &usize| percentage(*n, *total);
                let body = format!(
                    "\n\nAllele-specific paired-end sorting report\n{}\nRead pairs/singletons processed in total:\t\t{total}\n\tthereof were read pairs:\t\t\t{pairs}\n\tthereof were singletons:\t\t\t{singletons}\nReads were unassignable (not overlapping SNPs):\t\t{unassigned} ({}%)\n\tthereof were read pairs:\t{unassigned_pairs}\n\tthereof were singletons:\t{unassigned_singletons}\nReads were specific for genome 1:\t\t\t{genome1} ({}%)\n\tthereof were read pairs:\t{g1_pairs}\n\tthereof were singletons:\t{g1_singletons}\nReads were specific for genome 2:\t\t\t{genome2} ({}%)\n\tthereof were read pairs:\t{g2_pairs}\n\tthereof were singletons:\t{g2_singletons}\nReads contained conflicting SNP information:\t\t{conflicting} ({}%)\n\tthereof were read pairs:\t{conflicting_pairs}\n\tthereof were singletons:\t{conflicting_singletons}\n\n",
                    "=".repeat(41),
                    p(unassigned),
                    p(genome1),
                    p(genome2),
                    p(conflicting),
                );
                write!(err, "{body}")?;
                write!(outputs.report, "{body}")?;

                outputs.yaml("PE_total_reads", total);
                outputs.yaml("PE_total_pairs", pairs);
                outputs.yaml("PE_total_singletons", singletons);
                outputs.yaml("PE_unassignable", unassigned);
                outputs.yaml("PE_percent_unassignable", p(unassigned));
                outputs.yaml("PE_unassignable_pairs", unassigned_pairs);
                outputs.yaml("PE_unassignable_singletons", unassigned_singletons);
                outputs.yaml("PE_genome1", genome1);
                outputs.yaml("PE_percent_genome1", p(genome1));
                outputs.yaml("PE_genome1_pairs", g1_pairs);
                outputs.yaml("PE_genome1_singletons", g1_singletons);
                outputs.yaml("PE_genome2", genome2);
                outputs.yaml("PE_percent_genome2", p(genome2));
                outputs.yaml("PE_genome2_pairs", g2_pairs);
                outputs.yaml("PE_genome2_singletons", g2_singletons);
                outputs.yaml("PE_conflicting", conflicting);
                outputs.yaml("PE_percent_conflicting", p(conflicting));
                outputs.yaml("PE_conflicting_pairs", conflicting_pairs);
                outputs.yaml("PE_conflicting_singletons", conflicting_singletons);
            }

            Counts::HiC {
                total,
                unassigned,
                genome1,
                genome2,
                conflicting,
                g1_ua,
                ua_g1,
                g2_ua,
                ua_g2,
                g1_g2,
                g2_g1,
            } => {
                let p = |n: usize| percentage(n, *total);
                let g1_ua_total = g1_ua + ua_g1;
                let g2_ua_total = g2_ua + ua_g2;
                let g1_g2_total = g1_g2 + g2_g1;

                // The two copies differ by one word in the heading: stderr says "paired-end",
                // the report says "Hi-C paired-end".
                let body = |heading: &str| {
                    format!(
                        "\n\nAllele-specific {heading} sorting report\n{}\nRead pairs processed in total:\t\t\t\t{total}\nRead pairs were unassignable (UA/UA):\t\t\t{unassigned} ({}%)\nRead pairs were specific for genome 1 (G1/G1):\t\t{genome1} ({}%)\nRead pairs were specific for genome 2 (G2/G2):\t\t{genome2} ({}%)\nRead pairs were a mix of G1 and UA:\t\t\t{g1_ua_total} ({}%). Of these,\n\t\t\twere G1/UA: {g1_ua}\n\t\t\twere UA/G1: {ua_g1}\nRead pairs were a mix of G2 and UA:\t\t\t{g2_ua_total} ({}%). Of these,\n\t\t\twere G2/UA: {g2_ua}\n\t\t\twere UA/G2: {ua_g2}\nRead pairs were a mix of G1 and G2:\t\t\t{g1_g2_total} ({}%). Of these,\n\t\t\twere G1/G2: {g1_g2}\n\t\t\twere G2/G1: {g2_g1}\nRead pairs contained conflicting SNP information:\t{conflicting} ({}%)\n\n",
                        "=".repeat(41),
                        p(*unassigned),
                        p(*genome1),
                        p(*genome2),
                        p(g1_ua_total),
                        p(g2_ua_total),
                        p(g1_g2_total),
                        p(*conflicting),
                    )
                };
                write!(err, "{}", body("paired-end"))?;
                write!(outputs.report, "{}", body("Hi-C paired-end"))?;

                outputs.yaml("HiC_total_pairs", total);
                outputs.yaml("HiC_unassignable_UA_UA", unassigned);
                outputs.yaml("HiC_percent_unassignable_UA_UA", p(*unassigned));
                outputs.yaml("HiC_genome1_G1_G1", genome1);
                outputs.yaml("HiC_percent_genome1_G1_G1", p(*genome1));
                outputs.yaml("HiC_genome2_G2_G2", genome2);
                outputs.yaml("HiC_percent_genome2_G2_G2", p(*genome2));
                outputs.yaml("HiC_G1_UA_total", g1_ua_total);
                outputs.yaml("HiC_percent_G1_UA_total", p(g1_ua_total));
                outputs.yaml("HiC_G1_UA", g1_ua);
                outputs.yaml("HiC_UA_G1", ua_g1);
                outputs.yaml("HiC_G2_UA_total", g2_ua_total);
                outputs.yaml("HiC_percent_G2_UA_total", p(g2_ua_total));
                outputs.yaml("HiC_G2_UA", g2_ua);
                outputs.yaml("HiC_UA_G2", ua_g2);
                outputs.yaml("HiC_G1_G2_total", g1_g2_total);
                outputs.yaml("HiC_percent_G1_G2_total", p(g1_g2_total));
                outputs.yaml("HiC_G1_G2", g1_g2);
                outputs.yaml("HiC_G2_G1", g2_g1);
                outputs.yaml("HiC_conflicting", conflicting);
                outputs.yaml("HiC_percent_conflicting", p(*conflicting));
            }
        }
        Ok(())
    }
}
