//! Command line handling for `SNPsplit`, including the library-type autodetection.

use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::Result;

use crate::io::RecordReader;

/// The option names this tool registers, in the order the Perl's GetOptions lists them.
///
/// Resolution goes through `optmatch`, which reproduces Getopt::Long's case-insensitivity
/// and prefix matching, so every spelling an existing script might use still works.
/// Options this port adds. A tie with an option the Perl already had goes to the
/// original, so no existing abbreviation changes meaning.
const ADDED: &[&str] = &["parallel"];

pub const NAMES: &[&str] = &[
    "help",
    "man",
    "versions",
    "output_dir",
    "SNP_file",
    "no_sorting",
    "verbose",
    "samtools_path",
    "sam",
    "paired",
    "single_end",
    "hic",
    "conflicting",
    "weird",
    "singletons",
    "bisulfite",
    "skip_tag2sort",
    "parallel",
];

/// Raw options, before validation.
#[derive(Debug, Default)]
pub struct Options {
    pub snp_file: Option<PathBuf>,
    pub no_sort: bool,
    pub verbose: bool,
    pub samtools_path: Option<String>,
    pub sam: bool,
    pub paired: bool,
    pub single_end: bool,
    pub hic: bool,
    pub conflict: bool,
    pub output_dir: Option<String>,
    pub singletons: bool,
    pub bisulfite: bool,
    pub skip_tag2sort: bool,
    pub help: bool,
    pub version: bool,
    pub parallel: usize,
    pub files: Vec<PathBuf>,
}

/// The resolved configuration.
pub struct Config {
    pub snp_file: PathBuf,
    pub no_sort: bool,
    pub verbose: bool,
    pub samtools: String,
    pub samtools_given: Option<String>,
    pub bam: bool,
    pub paired: bool,
    pub hic: bool,
    pub conflict: bool,
    pub output_dir: String,
    pub parent_dir: String,
    pub singletons: bool,
    pub bisulfite: bool,
    pub skip_tag2sort: bool,
    pub parallel: usize,
    pub files: Vec<PathBuf>,
}

pub fn parse(args: &[String]) -> Result<Options, String> {
    let mut opts = Options {
        parallel: crate::io::default_parallel(),
        ..Options::default()
    };
    let mut i = 0;

    while i < args.len() {
        let arg = &args[i];
        let Some(body) = arg.strip_prefix("--").or_else(|| arg.strip_prefix("-")) else {
            opts.files.push(PathBuf::from(arg));
            i += 1;
            continue;
        };

        let (spelled, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v.to_string())),
            None => (body, None),
        };
        let name = crate::optmatch::resolve_with(spelled, NAMES, ADDED)?;

        let take_value = |i: &mut usize| -> Result<String, String> {
            if let Some(v) = inline.clone() {
                return Ok(v);
            }
            *i += 1;
            args.get(*i)
                .cloned()
                .ok_or_else(|| format!("Option {name} requires an argument"))
        };

        match name {
            "help" | "man" => opts.help = true,
            // Plural here, singular in tag2sort. See src/tool.rs.
            "versions" => opts.version = true,
            "output_dir" => opts.output_dir = Some(take_value(&mut i)?),
            "SNP_file" => opts.snp_file = Some(PathBuf::from(take_value(&mut i)?)),
            "no_sorting" => opts.no_sort = true,
            "verbose" => opts.verbose = true,
            "samtools_path" => opts.samtools_path = Some(take_value(&mut i)?),
            "sam" => opts.sam = true,
            "paired" => opts.paired = true,
            "single_end" => opts.single_end = true,
            "hic" => opts.hic = true,
            "conflicting" | "weird" => opts.conflict = true,
            "singletons" => opts.singletons = true,
            "bisulfite" => opts.bisulfite = true,
            "skip_tag2sort" => opts.skip_tag2sort = true,
            "parallel" => opts.parallel = take_value(&mut i)?.parse().unwrap_or(1),
            other => return Err(format!("Unknown option: {other}")),
        }

        i += 1;
    }

    Ok(opts)
}

/// Validate the options and work out the library type.
pub fn resolve(mut opts: Options) -> Result<Option<Config>> {
    if opts.files.is_empty() {
        // The warning alone is not actionable, which is why the Perl follows it with the
        // whole help page.
        eprintln!(
            "You need to provide one or more SAM/BAM files to start the allele-specific pipeline. Please respecify!"
        );
        print!("{}", crate::help::page(crate::Tool::Tag));
        anyhow::bail!("");
    }

    let Some(snp_file) = opts.snp_file.clone() else {
        anyhow::bail!(
            "You need to provide a text file detailing SNPs with '--SNP_file your.file'. Please respecify!\n"
        );
    };

    if opts.no_sort {
        for file in &opts.files {
            if !file.to_string_lossy().ends_with(".bam") {
                anyhow::bail!(
                    "The option --no_sorting requires you to supply a BAM file (ending in .bam), but at least one of the file names did not look like a BAM file: '{}'. Please respecify!\n\n",
                    file.display()
                );
            }
        }
    }

    if opts.single_end {
        if opts.paired {
            anyhow::bail!(
                "You cannot spcecify --paired as well as --single_end at the same time. Make your pick and try again!\n\n"
            );
        }
        opts.paired = false;
    }

    for file in &opts.files {
        if !file.exists() {
            anyhow::bail!(
                "Input file '{}' doesn't exist in the input folder. Please check filenames and try again!\n\n",
                file.display()
            );
        }
    }

    let samtools = resolve_samtools(opts.samtools_path.as_deref())?;

    if !snp_file.exists() {
        anyhow::bail!(
            "The SNP file '{}' does not exist. Please respecify!\n\n",
            snp_file.display()
        );
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

    let mut paired = opts.paired;
    let mut no_sort = opts.no_sort;

    if opts.hic {
        eprintln!(
            "Hi-C data specified. This assumes that the data has been processed with HiCUP, is paired-end and does not require positional sorting.\nSetting '--paired' and '--no_sort'...\n"
        );
        paired = true;
        no_sort = true;
    }

    let parent_dir = std::env::current_dir()?.to_string_lossy().to_string();
    let output_dir = match opts.output_dir.take() {
        None => String::new(),
        Some(dir) if dir.is_empty() => String::new(),
        Some(dir) => {
            let dir = if dir.ends_with('/') {
                dir
            } else {
                format!("{dir}/")
            };
            if !Path::new(&dir).is_dir() {
                std::fs::create_dir_all(&dir)
                    .map_err(|e| anyhow::anyhow!("Unable to create directory {dir} {e}\n"))?;
                eprintln!("Created output directory {dir}!\n");
            }
            let absolute = std::fs::canonicalize(&dir)?.to_string_lossy().to_string();
            let absolute = if absolute.ends_with('/') {
                absolute
            } else {
                format!("{absolute}/")
            };
            eprintln!("Output will be written into the directory: {absolute}");
            absolute
        }
    };

    let mut bisulfite = opts.bisulfite;

    if !opts.hic {
        let first = &opts.files[0];
        eprintln!(
            "Testing if input file '{}' looks like a Bisulfite-Seq file",
            first.display()
        );
        let (is_bismark, detected_paired) = check_for_bs(first, opts.paired, opts.single_end)?;
        paired = detected_paired;

        if is_bismark {
            bisulfite = true;
            if paired {
                eprintln!(
                    "File looks like a Bismark paired-end file. Setting '--bisulfite' and '--paired'..."
                );
                if positionally_sorted_is_fine(first)? {
                    eprintln!(
                        "File appears to be in the right format (Read1 and Read2 following each other), setting '--no_sort'..."
                    );
                    no_sort = true;
                } else {
                    eprintln!(
                        "File needs sorting by name first (Read1 and Read2 are required to follow each other for SNPsplit processing)..."
                    );
                }
            } else {
                eprintln!("File looks like a Bismark single-end file. Setting '--bisulfite'...");
            }
        }
    }

    Ok(Some(Config {
        snp_file,
        no_sort,
        verbose: opts.verbose,
        samtools,
        samtools_given: opts.samtools_path,
        bam: !opts.sam,
        paired,
        hic: opts.hic,
        conflict: opts.conflict,
        output_dir,
        parent_dir,
        singletons: opts.singletons,
        bisulfite,
        skip_tag2sort: opts.skip_tag2sort,
        parallel: opts.parallel,
        files: opts.files,
    }))
}

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
    path.split(':')
        .map(|dir| Path::new(dir).join("samtools"))
        .find(|candidate| candidate.exists())
        .map(|candidate| candidate.to_string_lossy().to_string())
}

/// Read the header's `@PG` records to decide whether this is Bismark output, and whether the
/// library is paired-end.
///
/// The aligner records its own command line there, and `-1`/`-2` in it is what says the run
/// was paired-end. This is why the fixture runner keeps the aligner's `@PG` record: it is an
/// input to the tool, not decoration.
pub fn check_for_bs(file: &Path, paired: bool, single_end: bool) -> Result<(bool, bool)> {
    let reader = RecordReader::open(file)?;
    let header = reader.header();

    let mut is_bismark = false;
    let mut detected = paired;

    for (id, program) in header.programs().as_ref().iter() {
        let id = String::from_utf8_lossy(id.as_ref());
        if id.starts_with("Bismark") {
            is_bismark = true;
        }
        if single_end {
            detected = false;
            continue;
        }
        if !(id.starts_with("Bismark") || id.starts_with("hisat2") || id.starts_with("bowtie2")) {
            continue;
        }
        if paired {
            continue;
        }

        let command = program
            .other_fields()
            .iter()
            .find(|(tag, _)| tag.as_ref() == b"CL")
            .map(|(_, value)| String::from_utf8_lossy(value.as_ref()).to_string())
            .unwrap_or_default();

        if command.contains(" -1 ") && command.contains(" -2 ") {
            eprintln!("Treating file(s) as paired-end data (as extracted from @PG line)\n");
            detected = true;
        } else {
            eprintln!("Treating file(s) as single-end data (as extracted from @PG line)\n");
            detected = false;
        }
        break;
    }

    Ok((is_bismark, detected))
}

/// Whether read 1 and read 2 already follow each other, so no name sort is needed.
pub fn positionally_sorted_is_fine(file: &Path) -> Result<bool> {
    eprintln!(
        "\nNow testing Bismark result file {} for positional sorting",
        file.display()
    );
    let reader = RecordReader::open(file)?;
    let header = reader.header().clone();

    // A header that declares coordinate sorting settles it without reading records.
    if let Some(hd) = header.header()
        && let Some(order) = hd
            .other_fields()
            .iter()
            .find(|(tag, _)| tag.as_ref() == b"SO")
        && String::from_utf8_lossy(order.1.as_ref()) == "coordinate"
    {
        eprintln!(
            "SAM header line indicates that the Bismark aligment file has been sorted by chromosomal positions. Paired-end files will be sorted by name before proceeding with the SNPsplit tagging process\n"
        );
        return Ok(false);
    }

    let mut previous: Option<Vec<u8>> = None;
    let mut count = 0usize;
    for record in reader {
        let record = record?;
        let name: Vec<u8> = record
            .name()
            .map(|n| {
                let bytes: &[u8] = n.as_ref();
                bytes.to_vec()
            })
            .unwrap_or_default();
        count += 1;
        if count > 100_000 {
            break;
        }
        match previous.take() {
            None => previous = Some(name),
            Some(first) => {
                if trimmed(&first) != trimmed(&name) {
                    eprintln!(
                        "The IDs of Read 1 ({}) and Read 2 ({}) are not the same. This might be a result of sorting the paired-end SAM/BAM files by chromosomal position. Paired-end files will be sorted by name before proceeding with the SNPsplit tagging process\n",
                        String::from_utf8_lossy(&first),
                        String::from_utf8_lossy(&name)
                    );
                    return Ok(false);
                }
            }
        }
    }

    Ok(true)
}

/// Older Bismark appended `/1` and `/2` to read names for readability.
fn trimmed(name: &[u8]) -> &[u8] {
    if name.ends_with(b"/1") || name.ends_with(b"/2") {
        &name[..name.len() - 2]
    } else {
        name
    }
}

/// Unused import guard: BufRead is needed by callers of this module's helpers.
#[allow(dead_code)]
fn _assert_buf_read<R: BufRead>(_: R) {}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_the_shape_the_fixtures_use() {
        let o = parse(&args(&[
            "--SNP_file",
            "snps.txt",
            "--single_end",
            "input.bam",
        ]))
        .unwrap();
        assert_eq!(o.snp_file.unwrap().to_str().unwrap(), "snps.txt");
        assert!(o.single_end);
        assert_eq!(o.files.len(), 1);
    }

    #[test]
    fn older_bismark_read_name_suffixes_are_ignored_when_pairing() {
        assert_eq!(trimmed(b"read1/1"), b"read1");
        assert_eq!(trimmed(b"read1/2"), b"read1");
        assert_eq!(trimmed(b"read1"), b"read1");
    }
}
