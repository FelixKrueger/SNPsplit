//! `tag2sort`: split allele-tagged alignments into per-genome files.

pub mod cli;
mod outputs;
mod report;

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use anyhow::Result;
use noodles_sam::alignment::RecordBuf;

use crate::io::RecordReader;
use crate::version;

use outputs::Outputs;

/// The allele call a tagged read carries in its `XX:Z:` tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Allele {
    /// Unassignable: the read overlapped no SNP.
    Unassigned,
    /// Genome 1 specific.
    Genome1,
    /// Genome 2 specific.
    Genome2,
    /// Conflicting SNP evidence within one read.
    Conflicting,
}

/// Read the `XX:Z:` tag SNPsplit wrote.
pub fn allele_of(record: &RecordBuf) -> Option<Allele> {
    let value = record
        .data()
        .get(&noodles_sam::alignment::record::data::field::Tag::new(
            b'X', b'X',
        ))?;
    // SNPsplit writes XX as a Z (string) field; anything else is not a tag we wrote.
    let noodles_sam::alignment::record_buf::data::field::Value::String(text) = value else {
        return None;
    };
    parse_allele(&String::from_utf8_lossy(text))
}

/// The four tags SNPsplit emits. Anything else is fatal, as in the Perl.
pub fn parse_allele(text: &str) -> Option<Allele> {
    match text {
        "UA" => Some(Allele::Unassigned),
        "G1" => Some(Allele::Genome1),
        "G2" => Some(Allele::Genome2),
        "CF" => Some(Allele::Conflicting),
        _ => None,
    }
}

/// Entry point.
pub fn run(args: &[String]) -> ExitCode {
    match execute(args) {
        Ok(code) => code,
        Err(e) => {
            eprint!("{e}");
            ExitCode::FAILURE
        }
    }
}

fn execute(args: &[String]) -> Result<ExitCode> {
    let mut opts = match cli::parse(args) {
        Ok(o) => o,
        Err(message) => {
            eprintln!("{message}");
            anyhow::bail!("Please respecify command line options\n");
        }
    };

    if opts.version {
        print!("{}", version::banner(crate::Tool::Sort));
        return Ok(ExitCode::SUCCESS);
    }

    if opts.files.is_empty() {
        eprintln!(
            "You need to provide one or more allele-tagged SNPsplit files to start sorting them allele-specifically. Please respecify!\n"
        );
        return Ok(ExitCode::SUCCESS);
    }

    // An empty string is meaningful: SNPsplit passes one when it has no directory to add, and
    // turning it into "/" would send every output to the filesystem root.
    let output_dir = match opts.output_dir.take() {
        Some(dir) if dir.is_empty() => String::new(),
        Some(dir) if dir.ends_with('/') => dir,
        Some(dir) => format!("{dir}/"),
        None => String::new(),
    };

    for file in &opts.files {
        let path = format!("{output_dir}{}", file.display());
        if !Path::new(&path).exists() {
            anyhow::bail!(
                "Input file '{path}' doesn't exist in the input folder. Please check filenames and try again!\n\n"
            );
        }
        let name = file.to_string_lossy();
        if !(name.ends_with(".bam") || name.ends_with(".sam")) {
            anyhow::bail!(
                "Supplied file needs to be a SNPsplit tagged BAM or SAM file (ending in .bam or .sam). Please respecify!\n"
            );
        }
    }

    if opts.singletons {
        if opts.hic {
            anyhow::bail!(
                "The option --singletons can't be used in Hi-C mode since this expects paired-end data. Please respecify!\n\n"
            );
        }
        if !opts.paired {
            eprintln!(
                "The option --singletons only makes sense for paired-end files. Simply ignoring it in single-end mode..."
            );
            opts.singletons = false;
        }
    }

    let samtools = resolve_samtools(opts.samtools_path.as_deref())?;
    let parent_dir = std::env::current_dir()?.to_string_lossy().to_string();

    for file in &opts.files.clone() {
        sort_one(&opts, file, &output_dir, &parent_dir, &samtools)?;
    }

    Ok(ExitCode::SUCCESS)
}

/// Resolve `--samtools_path`, which no longer selects a reader but is still validated and
/// still reported.
///
/// The Rust build reads BAM itself, so this path is not run. It stays accepted so no existing
/// command line breaks, and it is still checked so a wrong one is still reported rather than
/// silently ignored. The reported value is masked by the fixture runner, so what it says
/// matters only to a human reading a log.
fn resolve_samtools(requested: Option<&str>) -> Result<String> {
    let Some(path) = requested else {
        return Ok(which_samtools().unwrap_or_else(|| "<not-required>".to_string()));
    };

    let candidate = if path.ends_with("samtools") {
        path.to_string()
    } else if path.ends_with('/') {
        format!("{path}samtools")
    } else {
        format!("{path}/samtools")
    };

    if Path::new(&candidate).exists() {
        Ok(candidate)
    } else {
        anyhow::bail!(
            "Could not find an installation of Samtools at the location {candidate}. Please respecify\n"
        );
    }
}

fn which_samtools() -> Option<String> {
    let path = std::env::var("PATH").ok()?;
    for dir in path.split(':') {
        let candidate = Path::new(dir).join("samtools");
        if candidate.exists() {
            return Some(candidate.to_string_lossy().to_string());
        }
    }
    None
}

fn sort_one(
    opts: &cli::Options,
    file: &Path,
    output_dir: &str,
    parent_dir: &str,
    samtools: &str,
) -> Result<()> {
    let mut err = std::io::stderr();
    let name = file.to_string_lossy().to_string();

    writeln!(err, "\nSummary of parameters for SNPsplit-sort:")?;
    writeln!(err, "{}", "=".repeat(40))?;
    writeln!(err, "SNPsplit tagged infile:\t\t{name}")?;
    writeln!(err, "Output directory:\t\t>{output_dir}<")?;
    writeln!(err, "Parent directory:\t\t>{parent_dir}<")?;
    writeln!(err, "Samtools path:\t\t\t{samtools}")?;
    writeln!(
        err,
        "Output format:\t\t\t{}",
        if opts.sam { "SAM" } else { "BAM (default)" }
    )?;
    if opts.hic {
        writeln!(err, "Input format:\t\t\tHi-C (by definition paired-end)")?;
    } else if opts.paired {
        if opts.singletons {
            writeln!(
                err,
                "Input format:\t\t\tPaired-end (Singleton alignments will written to extra files)"
            )?;
        } else {
            writeln!(err, "Input format:\t\t\tPaired-end")?;
        }
    } else {
        writeln!(err, "Input format:\t\t\tSingle-End")?;
    }
    writeln!(err, "\n")?;

    let input = format!("{output_dir}{name}");
    let reader = RecordReader::open(Path::new(&input))?;
    let mut header = reader.header().clone();
    crate::io::add_pg_line(
        &mut header,
        "SNPsplit",
        version::SUITE_VERSION,
        &crate::io::command_line(),
    )?;

    // Single-end announces the file before opening its outputs; the two paired modes do not.
    if !opts.hic && !opts.paired {
        writeln!(err, "Now processing input file <<< {name} >>>\n")?;
    }

    let mut outputs = Outputs::open(opts, &name, output_dir, &header, &mut err)?;
    // First YAML entry, recorded in the main flow before any processing.
    outputs.yaml("tagged_infile", &name);

    let counts = if opts.hic {
        report::process_hic(reader, &header, &mut outputs, &mut err)?
    } else if opts.paired {
        report::process_paired_end(reader, &header, &mut outputs, opts.singletons, &mut err)?
    } else {
        report::process_single_end(reader, &header, &mut outputs, opts.verbose, &mut err)?
    };

    // The report body is written by the processing step itself, before the run announces it
    // has finished.
    counts.write(&mut outputs, &mut err)?;
    outputs.finish()?;
    writeln!(err, "Sorting finished successfully\n")?;
    outputs.close_report()?;
    Ok(())
}

/// Percentages the Perl renders with `sprintf "%.2f"`, or `N/A` when nothing was counted.
pub fn percentage(part: usize, total: usize) -> String {
    if total == 0 {
        "N/A".to_string()
    } else {
        format!("{:.2}", part as f64 * 100.0 / total as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_the_four_snpsplit_tags_are_alleles() {
        assert_eq!(parse_allele("UA"), Some(Allele::Unassigned));
        assert_eq!(parse_allele("G1"), Some(Allele::Genome1));
        assert_eq!(parse_allele("G2"), Some(Allele::Genome2));
        assert_eq!(parse_allele("CF"), Some(Allele::Conflicting));
        assert_eq!(parse_allele("G3"), None);
        assert_eq!(parse_allele(""), None);
    }

    /// Nothing processed reports N/A rather than a division by zero, and the fixtures for
    /// header-only input diff exactly that.
    #[test]
    fn percentages_of_nothing_are_not_a_number_but_a_word() {
        assert_eq!(percentage(0, 0), "N/A");
        assert_eq!(percentage(1, 3), "33.33");
        assert_eq!(percentage(2, 3), "66.67");
        assert_eq!(percentage(1, 1), "100.00");
    }
}
