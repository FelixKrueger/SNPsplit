//! `SNPsplit_genome_preparation`: build N-masked and full-sequence genomes from a VCF.

pub mod cli;
pub mod download;
pub mod dual;
pub mod genome;
pub mod vcf;

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result};

use crate::version;

/// Entry point. Errors are printed the way Perl's `die` prints them and become a non-zero
/// exit, which is all the fixture suite asserts about the status.
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
    let opts = match cli::parse(args) {
        Ok(o) => o,
        Err(message) => {
            eprintln!("{message}");
            anyhow::bail!("Please respecify command line options\n");
        }
    };

    if opts.help {
        print!("{}", crate::help::page(crate::Tool::Prepare));
        return Ok(ExitCode::from(crate::help::exit_status(
            crate::Tool::Prepare,
        )));
    }

    if opts.version {
        print!("{}", version::banner(crate::Tool::Prepare));
        return Ok(ExitCode::SUCCESS);
    }

    let Some(config) = resolve(opts)? else {
        // --list_strains printed its list and stopped.
        return Ok(ExitCode::SUCCESS);
    };

    prepare(&config)?;
    Ok(ExitCode::SUCCESS)
}

/// Validate the options and look up the strain columns, mirroring `process_commandline`.
///
/// `Ok(None)` means the run is over: `--list_strains` prints and exits successfully.
fn resolve(mut opts: cli::Options) -> Result<Option<cli::Config>> {
    let mut v7 = opts.v7;

    // Before anything is validated, because what it fetches is what the validation then
    // looks at. Nothing below contacts the network when the flag is absent.
    if opts.download {
        fetch_inputs(&mut opts, v7)?;
    }

    if let Some(vcf) = &opts.vcf_file {
        if !vcf.exists() {
            anyhow::bail!(
                "Input VCF file '{}' doesn't exist in the folder. Please check filenames and try again!\n\n",
                vcf.display()
            );
        }
        if cli::is_v7_release(vcf) {
            eprintln!(
                "Setting option --v7_VCF as the supplied file is called 'mgp_REL2005_snps_indels.vcf.gz'\n"
            );
            v7 = true;
        }
        if v7 {
            if cli::basename(vcf) == "mgp_REL2005_snps_indels.vcf.gz" {
                eprintln!(
                    "Using v7 MGP file 'mgp_REL2005_snps_indels.vcf.gz'. This file can be obtained from:"
                );
                eprintln!(
                    "ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz\n"
                );
                eprintln!(
                    "PLEASE NOTE:\n============\nThis file is currently not marked as Current_SNPs, so consider this approach experimental for the time being (as of 15 03 2021)\n"
                );
            } else {
                eprintln!(
                    "Version 7 input file selected. Input VCF file '{}' doesn't appear to be file 'mgp_REL2005_snps_indels.vcf.gz'. If something goes wrong, I won't accept responsibility....!\n",
                    vcf.display()
                );
            }
        }
    } else if !opts.skip_filtering {
        anyhow::bail!(
            "\nYou need to provide a VCF file detailing SNPs positions with '--vcf_file your.file' (e.g.: --vcf mgp.v5.merged.snps_all.dbSNP142.vcf.gz). Please respecify!\n\n"
        );
    }

    let mut dual_hybrid = opts.dual_hybrid;
    let mut full_sequence = opts.full_sequence;

    if opts.strain2.is_some() && !dual_hybrid {
        eprintln!("Strain 2 specified, setting option '--dual_hybrid'");
        dual_hybrid = true;
    }

    if dual_hybrid {
        eprintln!("Dual Hybrid strain selected");
        full_sequence = true;

        let Some(strain2) = &opts.strain2 else {
            anyhow::bail!(
                "\n'--dual_hybrid' needs a second strain. Please respecify using '--strain2 NAME'\n\n"
            );
        };
        let Some(strain) = &opts.strain else {
            anyhow::bail!(
                "\n'--dual_hybrid' needs a first strain. Please respecify using '--strain NAME'\n\n"
            );
        };
        if strain == strain2 {
            anyhow::bail!(
                "Strain 1 [{strain}] and Strain 2 [{strain2}] must be different from each other. Please respecify!\n\n"
            );
        }
    }

    if opts.skip_filtering && opts.strain.is_none() {
        anyhow::bail!(
            "\nYou need to name the strain whose SNP folder should be used, with '--strain NAME'\n\n"
        );
    }

    let mut strain_index = None;
    let mut strain2_index = None;

    if !opts.skip_filtering {
        let vcf = opts.vcf_file.as_ref().expect("checked above");
        let strains = vcf::detect_strains(vcf)?;

        if opts.list_strains {
            eprintln!("\nAvailable genomes to choose from are:");
            eprintln!("{}", "=".repeat(37));
            for name in strains.keys() {
                println!("{name}");
            }
            eprintln!("{}", "=".repeat(37));
            eprintln!("\nPlease choose a strain using '--strain NAME' to continue.\n");
            return Ok(None);
        }

        match &opts.strain {
            Some(strain) => match strains.get(strain) {
                Some(index) => {
                    strain_index = Some(*index);
                    eprintln!("Strain defined as '{strain}' (strain index: {index})");
                }
                None => {
                    eprintln!(
                        "Strain name specified [{strain}] does not match any of the available strain names!"
                    );
                    list_available(&strains);
                    anyhow::bail!(
                        "\nPlease double check the name and try again (using '--strain NAME')\n\n"
                    );
                }
            },
            None => {
                eprintln!("No strain specified!");
                list_available(&strains);
                anyhow::bail!(
                    "\nPlease choose one of the available strains using '--strain NAME' and try again\n\n"
                );
            }
        }

        if dual_hybrid {
            match &opts.strain2 {
                Some(strain2) => match strains.get(strain2) {
                    Some(index) => {
                        strain2_index = Some(*index);
                        eprintln!("Strain2 defined as '{strain2}' (strain2 index: {index})");
                    }
                    None => {
                        eprintln!(
                            "Strain2 name specified [{strain2}] does not match any of the available strain names!"
                        );
                        list_available(&strains);
                        anyhow::bail!(
                            "\nPlease double check the name and try again (using '--strain2 NAME')\n\n"
                        );
                    }
                },
                None => {
                    eprintln!("No strain 2 specified!");
                    list_available(&strains);
                    anyhow::bail!(
                        "\nPlease choose one of the available strains using '--strain2 NAME' and try again\n\n"
                    );
                }
            }
        }
    }

    let Some(genome_folder) = opts.genome_folder else {
        anyhow::bail!(
            "Reference genome folder was not specified! Please use --reference_genome </genome/folder/>\n"
        );
    };

    let genome_folder = if genome_folder.ends_with('/') {
        genome_folder
    } else {
        format!("{genome_folder}/")
    };

    let absolute = match std::fs::canonicalize(&genome_folder) {
        Ok(p) => {
            let mut s = p.to_string_lossy().to_string();
            if !s.ends_with('/') {
                s.push('/');
            }
            s
        }
        Err(e) => {
            anyhow::bail!(
                "Failed to move to genome folder > {genome_folder} <: {e}\n\nSNPsplit_genome_preparation --help for more details\n\n"
            );
        }
    };
    eprintln!(
        "Reference genome folder provided is {genome_folder}\t(absolute path is '{absolute})'\n"
    );

    // --no_nmasking swaps the default rather than merely clearing it: a run has to produce
    // one genome or the other.
    let nmasking = !opts.no_nmasking;
    if opts.no_nmasking {
        full_sequence = true;
    }

    Ok(Some(cli::Config {
        vcf_file: opts.vcf_file,
        strain: opts.strain.unwrap_or_default(),
        strain2: opts.strain2,
        strain_index,
        strain2_index,
        genome_folder: absolute,
        skip_filtering: opts.skip_filtering,
        nmasking,
        full_sequence,
        dual_hybrid,
        genome_build: opts.genome_build.unwrap_or_else(|| "GRCm39".to_string()),
        v7,
        parallel: opts.parallel,
    }))
}

/// Fetch whatever of the two inputs is missing.
///
/// A path the user named always wins: if `--vcf_file` or `--reference_genome` points at
/// something that exists, it is used as given and nothing is downloaded over it. This is the
/// one place in the port that reaches the network, and only with `--download`.
fn fetch_inputs(opts: &mut cli::Options, v7: bool) -> Result<()> {
    let mut err = std::io::stderr();
    let dir = PathBuf::from(
        opts.download_dir
            .clone()
            .unwrap_or_else(|| "SNPsplit_references".to_string()),
    );
    let build = opts
        .genome_build
        .clone()
        .unwrap_or_else(|| "GRCm39".to_string());
    let release = opts.ensembl_release.unwrap_or(DEFAULT_ENSEMBL_RELEASE);
    let sources = download::sources(v7, &build, release);

    std::fs::create_dir_all(&dir)
        .with_context(|| format!("Failed to create download directory {}", dir.display()))?;
    writeln!(err, "Downloading reference data into {}\n", dir.display())?;

    let mut fetched = Vec::new();

    let wanted_vcf = opts
        .vcf_file
        .clone()
        .unwrap_or_else(|| dir.join(&sources.vcf_name));
    if wanted_vcf.exists() {
        writeln!(
            err,
            "Using the VCF already present at '{}'",
            wanted_vcf.display()
        )?;
    } else {
        let bytes = download::fetch(&sources.vcf_url, &wanted_vcf, &mut err)?;
        download::verify_vcf(&wanted_vcf)?;
        writeln!(err, "Verified '{}'", wanted_vcf.display())?;
        fetched.push(download::Fetched {
            url: sources.vcf_url.clone(),
            path: wanted_vcf.clone(),
            bytes,
            sha256: download::sha256(&wanted_vcf)?,
        });
    }
    opts.vcf_file = Some(wanted_vcf);

    let wanted_genome = PathBuf::from(
        opts.genome_folder
            .clone()
            .unwrap_or_else(|| dir.join(&build).to_string_lossy().to_string()),
    );
    if wanted_genome.is_dir() && has_fasta(&wanted_genome) {
        writeln!(
            err,
            "Using the reference genome already present at '{}'",
            wanted_genome.display()
        )?;
    } else {
        std::fs::create_dir_all(&wanted_genome)?;
        for chr in download::mouse_chromosomes() {
            let name = format!("Mus_musculus.{build}.dna.chromosome.{chr}.fa.gz");
            let target = wanted_genome.join(format!("{chr}.fa"));
            if target.exists() {
                continue;
            }
            let url = format!("{}/{name}", sources.genome_base);
            let compressed = wanted_genome.join(&name);
            let bytes = download::fetch(&url, &compressed, &mut err)?;
            decompress(&compressed, &target)?;
            fetched.push(download::Fetched {
                url,
                path: target,
                bytes,
                sha256: download::sha256(&compressed)?,
            });
            std::fs::remove_file(&compressed)?;
        }
    }
    opts.genome_folder = Some(wanted_genome.to_string_lossy().to_string());

    if !fetched.is_empty() {
        download::record(&dir, &fetched)?;
        writeln!(
            err,
            "\nRecorded {} download(s) in {}",
            fetched.len(),
            dir.join("manifest.txt").display()
        )?;
    }
    writeln!(err)?;
    Ok(())
}

/// The Ensembl release the genome is taken from when none is given.
const DEFAULT_ENSEMBL_RELEASE: u32 = 115;

fn has_fasta(dir: &Path) -> bool {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return false;
    };
    entries.filter_map(|e| e.ok()).any(|e| {
        let name = e.file_name().to_string_lossy().to_string();
        name.ends_with(".fa") || name.ends_with(".fasta")
    })
}

/// Decompress a downloaded chromosome, since the genome folder holds plain FastA.
fn decompress(source: &Path, target: &Path) -> Result<()> {
    let input = std::fs::File::open(source)?;
    let mut reader = flate2::read::MultiGzDecoder::new(input);
    let mut output = std::fs::File::create(target)
        .with_context(|| format!("Failed to write to {}", target.display()))?;
    std::io::copy(&mut reader, &mut output)
        .with_context(|| format!("'{}' is truncated or corrupt", source.display()))?;
    Ok(())
}

fn list_available(strains: &BTreeMap<String, usize>) {
    eprintln!("\nAvailable genomes to choose from are:");
    eprintln!("{}", "=".repeat(37));
    for name in strains.keys() {
        println!("{name}");
    }
    eprintln!("{}", "=".repeat(37));
}

/// Which of the three genome-building passes is running.
///
/// They differ in more than the strain name: the leading blank lines of the summary, the
/// label it carries, whether the per-chromosome decisions are announced, and, in the third
/// pass, the separators the report omits.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Pass {
    First,
    Second,
    DualHybrid,
}

/// The run proper: filter, read the genome, write the modified chromosomes.
fn prepare(config: &cli::Config) -> Result<()> {
    let parent = std::env::current_dir()?;
    let mut out = std::io::stdout();
    let mut err = std::io::stderr();

    summarise(config, &mut err)?;

    let chroms: Vec<String> = if config.skip_filtering {
        let mut v: Vec<String> = (1..=19).map(|n| n.to_string()).collect();
        v.extend(["X", "Y", "MT"].iter().map(|s| s.to_string()));
        v
    } else {
        vcf::detect_chroms(config.vcf_file.as_ref().expect("checked in resolve"))?
    };
    let chrom_set: BTreeSet<String> = chroms.iter().cloned().collect();

    if config.skip_filtering {
        println!("Using the following chromosomes (HARDCODED IN!!!):");
    } else {
        println!(
            "Using the following chromosomes (detected from VCF file >>{}<<):",
            config.vcf_file.as_ref().expect("checked").display()
        );
    }
    println!("{}\n", chroms.join("\t"));

    // Shared across both filtering passes: the dual hybrid comparison needs to tell a
    // position that was absent from one that was present but low confidence.
    let mut blacklist: BTreeMap<String, vcf::Blacklisted> = BTreeMap::new();

    if config.skip_filtering {
        writeln!(
            err,
            "Skipped reading the VCF file and filtering SNPs again (specified by user)\n"
        )?;
    } else {
        vcf::filter_snps(
            config.vcf_file.as_ref().expect("checked"),
            &config.strain,
            config.strain_index.expect("checked"),
            vcf::StrainIdentity::First,
            &chroms,
            &config.genome_build,
            config.v7,
            config.dual_hybrid,
            &mut blacklist,
            &mut out,
            &mut err,
        )?;
        writeln!(
            err,
            "Finished filtering and writing out SNPs for strain {}\n",
            config.strain
        )?;
    }

    let genome =
        genome::read_genome_into_memory(Path::new(&config.genome_folder), &mut out, &mut err)?;

    check_chromosome_names_agree(&genome, &chrom_set, &mut err)?;

    build_genome(
        config,
        &parent,
        &genome,
        &chrom_set,
        Pass::First,
        &config.strain.clone(),
        None,
        &mut err,
    )?;

    if config.dual_hybrid {
        let strain2 = config.strain2.clone().expect("checked in resolve");
        writeln!(err, "Now starting to work on strain 2 [{strain2}]")?;

        if config.skip_filtering {
            writeln!(
                err,
                "Skipped reading the VCF file and filtering SNPs again for strain 2 (specified by user)\n"
            )?;
        } else {
            vcf::filter_snps(
                config.vcf_file.as_ref().expect("checked"),
                &strain2,
                config.strain2_index.expect("checked"),
                vcf::StrainIdentity::Second,
                &chroms,
                &config.genome_build,
                config.v7,
                config.dual_hybrid,
                &mut blacklist,
                &mut out,
                &mut err,
            )?;
            writeln!(
                err,
                "Finished filtering and writing out SNPs for strain 2 [{strain2}]\n"
            )?;
        }

        build_genome(
            config,
            &parent,
            &genome,
            &chrom_set,
            Pass::Second,
            &strain2,
            None,
            &mut err,
        )?;

        build_dual_hybrid(
            config, &parent, &chrom_set, &strain2, &blacklist, &mut out, &mut err,
        )?;
    }

    writeln!(
        err,
        "All done. Genome(s) are now ready to be indexed with your favourite aligner!\nFYI, aligners shown to work with SNPsplit are Bowtie2, STAR, HISAT2, HiCUP and Bismark (STAR and Hisat2 require disabling soft-clipping, please check the SNPsplit manual for details)\n"
    )?;

    Ok(())
}

/// The third genome: strain 1's full sequence as the reference, strain 2's bases as the SNPs.
fn build_dual_hybrid(
    config: &cli::Config,
    parent: &Path,
    chroms: &BTreeSet<String>,
    strain2: &str,
    blacklist: &BTreeMap<String, vcf::Blacklisted>,
    out: &mut impl Write,
    err: &mut impl Write,
) -> Result<()> {
    let strain = &config.strain;
    let build = &config.genome_build;

    let report_name = format!("{strain}_{strain2}_dual_hybrid.genome_preparation_report.txt");
    let mut report = std::io::BufWriter::new(std::fs::File::create(&report_name)?);

    // Named relative to the working directory, as the Perl does: the name is echoed into the
    // log and an absolute path there would differ from the recorded output.
    let snp_file_strain = PathBuf::from(format!("all_SNPs_{strain}_{build}.txt.gz"));
    let snp_file_strain2 = PathBuf::from(format!("all_SNPs_{strain2}_{build}.txt.gz"));

    let (annotation, _counts) = dual::determine_snps_between_strains(
        strain,
        strain2,
        build,
        &snp_file_strain,
        &snp_file_strain2,
        blacklist,
        &mut report,
        err,
    )?;
    writeln!(err, "done...")?;

    let line = format!(
        "Changing the genomic reference sequence to the full sequence of strain {strain}\n"
    );
    writeln!(err, "{line}")?;
    writeln!(report, "{line}")?;

    // The trailing slash is part of the logged path: the Perl appends one to every genome
    // folder before reporting it.
    let full_folder = parent.join(format!("{strain}_full_sequence/"));
    let genome = genome::read_genome_into_memory(&full_folder, out, err)?;

    let line = format!(
        "Reading and storing all new SNPs with Ref/SNP: {strain}/{strain2} from '{annotation}'"
    );
    writeln!(err, "{line}")?;
    writeln!(report, "{line}")?;

    let snps = dual::read_new_snp_annotation(Path::new(&annotation), strain, strain2, err)?;

    // build_genome reopens this file to append its own half, so the buffered writes have to
    // land first. A BufWriter still holding bytes would flush at its own offset afterwards
    // and overwrite what the second half wrote.
    report.flush()?;
    drop(report);

    let label = format!("{strain}_{strain2}_dual_hybrid.based_on_{build}");
    build_genome(
        config,
        parent,
        &genome,
        chroms,
        Pass::DualHybrid,
        &label,
        Some(&snps),
        err,
    )?;

    // The Perl writes these two summaries to the report without the separators its stderr
    // copy carries, so the file says "strainstrain 2 [AB]" where the terminal says
    // "strain/strain 2 [A/B]". Reproduced; raised upstream separately.
    Ok(())
}

/// Build one genome: apply the SNPs each chromosome has, write the rest through unchanged.
#[allow(clippy::too_many_arguments)]
fn build_genome(
    config: &cli::Config,
    parent: &Path,
    genome: &genome::Genome,
    chroms: &BTreeSet<String>,
    pass: Pass,
    label: &str,
    preloaded: Option<&BTreeMap<String, Vec<genome::Snp>>>,
    err: &mut impl Write,
) -> Result<()> {
    let report_name = match pass {
        Pass::First => format!("{}_genome_preparation_report.txt", config.strain),
        Pass::Second => format!("{label}_genome_preparation_report.txt"),
        // The dual pass appends to the report the caller already opened, which is reopened
        // here in append mode so both halves land in one file.
        Pass::DualHybrid => format!(
            "{}_{}_dual_hybrid.genome_preparation_report.txt",
            config.strain,
            config.strain2.as_deref().unwrap_or("")
        ),
    };
    let mut report = std::io::BufWriter::new(
        std::fs::OpenOptions::new()
            .create(true)
            .append(pass == Pass::DualHybrid)
            .write(true)
            .truncate(pass != Pass::DualHybrid)
            .open(&report_name)?,
    );

    let mut new_n_total = 0usize;
    let mut new_snp_total = 0usize;

    // Chromosomes are independent: each reads its own SNP track and writes its own file. The
    // work runs in parallel and the output is emitted in genome order afterwards, so the
    // logs, the report and the totals are the same whatever the worker count. An abort
    // surfaces in genome order too, so which chromosome stops the run does not depend on
    // which worker reached it first.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(crate::io::worker_count(config.parallel).get())
        .build()?;

    let outcomes: Vec<ChromosomeOutcome> = pool.install(|| {
        use rayon::prelude::*;
        genome
            .par_iter()
            .map(|(chr, sequence)| {
                one_chromosome(
                    config, parent, genome, chroms, pass, label, preloaded, chr, sequence,
                )
            })
            .collect()
    });

    for outcome in outcomes {
        err.write_all(&outcome.log)?;
        if let Some(e) = outcome.error {
            return Err(e);
        }
        report.write_all(&outcome.report)?;
        for (path, chr, sequence) in &outcome.writes {
            genome::write_planned(path, chr, sequence)?;
        }
        new_n_total += outcome.new_n;
        new_snp_total += outcome.new_snp;
    }

    write_summary(
        config,
        pass,
        label,
        new_n_total,
        new_snp_total,
        &mut report,
        err,
    )?;
    report.flush()?;
    Ok(())
}

/// What one chromosome produced: its log, its report lines, and its counters.
struct ChromosomeOutcome {
    log: Vec<u8>,
    report: Vec<u8>,
    new_n: usize,
    new_snp: usize,
    /// Files this chromosome wants written, performed by the caller in genome order so an
    /// earlier abort leaves none of them behind.
    writes: Vec<(std::path::PathBuf, String, Vec<u8>)>,
    /// An abort, carried rather than returned so it surfaces during the ordered replay. A
    /// worker that returned it directly would lose the lines it had already logged, and would
    /// report the failure before earlier chromosomes had reported their success.
    error: Option<anyhow::Error>,
}

/// Process one chromosome, buffering everything it would have printed.
///
/// Buffering is what makes the parallelism invisible: the caller replays these in genome
/// order, so the log and the report read exactly as they would from a serial run.
#[allow(clippy::too_many_arguments)]
fn one_chromosome(
    config: &cli::Config,
    parent: &Path,
    genome: &genome::Genome,
    chroms: &BTreeSet<String>,
    pass: Pass,
    label: &str,
    preloaded: Option<&BTreeMap<String, Vec<genome::Snp>>>,
    chr: &str,
    sequence: &[u8],
) -> ChromosomeOutcome {
    let mut log: Vec<u8> = Vec::new();
    let mut report: Vec<u8> = Vec::new();
    let mut writes: Vec<(std::path::PathBuf, String, Vec<u8>)> = Vec::new();

    let result = chromosome_body(
        config,
        parent,
        genome,
        chroms,
        pass,
        label,
        preloaded,
        chr,
        sequence,
        &mut log,
        &mut report,
        &mut writes,
    );

    match result {
        Ok((new_n, new_snp)) => ChromosomeOutcome {
            log,
            report,
            new_n,
            new_snp,
            writes,
            error: None,
        },
        Err(e) => ChromosomeOutcome {
            log,
            report,
            new_n: 0,
            new_snp: 0,
            // An aborting chromosome leaves nothing behind.
            writes: Vec::new(),
            error: Some(e),
        },
    }
}

#[allow(clippy::too_many_arguments)]
fn chromosome_body(
    config: &cli::Config,
    parent: &Path,
    genome: &genome::Genome,
    chroms: &BTreeSet<String>,
    pass: Pass,
    label: &str,
    preloaded: Option<&BTreeMap<String, Vec<genome::Snp>>>,
    chr: &str,
    sequence: &[u8],
    log: &mut Vec<u8>,
    report: &mut Vec<u8>,
    writes: &mut Vec<(std::path::PathBuf, String, Vec<u8>)>,
) -> Result<(usize, usize)> {
    if !chroms.contains(chr) {
        let err = &mut *log;
        if pass == Pass::DualHybrid {
            writeln!(
                err,
                "Got no SNP information for chromosome {chr}. Printing sequence only..."
            )?;
        }
        if config.nmasking {
            let path = genome::plan_chromosome(parent, chr, true, label, err)?;
            writes.push((path, chr.to_string(), sequence.to_vec()));
        }
        if config.full_sequence {
            let path = genome::plan_chromosome(parent, chr, false, label, err)?;
            writes.push((path, chr.to_string(), sequence.to_vec()));
        }
        return Ok((0, 0));
    }

    {
        let err = &mut *log;
        match pass {
            Pass::DualHybrid => writeln!(
                err,
                "Got SNP information for chromosome {chr}. Creating modified chromosome"
            )?,
            _ => writeln!(err, "Processing chromosome {chr} (for strain {label})")?,
        }

        if pass == Pass::DualHybrid {
            writeln!(
                err,
                "Processing chr{chr} (to create new genome for {}/{})",
                config.strain,
                config.strain2.as_deref().unwrap_or("")
            )?;
        }
    }

    // The Perl guards this loop with `unless ($chromosomes{$chr})`, and an empty string is
    // false in Perl, so a chromosome that exists but carries no sequence is reported as a
    // name mismatch. Reproduced; the misleading diagnostic is raised upstream separately.
    if sequence.is_empty() && pass != Pass::DualHybrid {
        let err = &mut *log;
        writeln!(
            err,
            "\nThe chromosome name given in the VCF file was '{chr}' and was not found in the reference genome.\nA rather common mistake might be that the VCF file was downloaded from Ensembl (who use chromosome names such as 1, 2, X, MT)\nbut the genome from UCSC (who use chromosome names such as chr1, chr2, chrX, chrM)"
        )?;
        writeln!(
            err,
            "The chromosome names in the reference genome folder were:"
        )?;
        for name in genome.keys() {
            writeln!(err, "{name}")?;
        }
        anyhow::bail!(
            "[FATAL ERROR] Please ensure that the same version of the genome is used for both VCF annotations and reference genome (FastA files). Exiting...\n\n"
        );
    }

    let err = &mut *log;
    let owned;
    let snps: &[genome::Snp] = match preloaded {
        Some(map) => map.get(chr).map(Vec::as_slice).unwrap_or(&[]),
        None => {
            owned = genome::read_snps(parent, chr, label, err)?;
            if owned.is_empty() {
                // Assigning an empty list to a list that is already empty, announced.
                writeln!(err, "Clearing SNP array...")?;
            }
            &owned
        }
    };

    let (counts, masked, full) =
        genome::apply_snps(sequence, snps, config.nmasking, config.full_sequence);

    if let Some(masked) = masked {
        let path = genome::plan_chromosome(parent, chr, true, label, err)?;
        writes.push((path, chr.to_string(), masked));
    }
    if let Some(full) = full {
        let path = genome::plan_chromosome(parent, chr, false, label, err)?;
        writes.push((path, chr.to_string(), full));
    }

    writeln!(err, "{} SNPs total for chromosome {chr}", counts.total)?;
    if counts.already_carried > 0 {
        let line = format!(
            "{} positions on chromosome {chr} already carried the SNP base and were left alone",
            counts.already_carried
        );
        writeln!(err, "{line}")?;
        writeln!(report, "{line}")?;
    }
    if counts.mismatched > 0 {
        let line = format!(
            "{} positions on chromosome {chr} were skipped because the reference base did not match the annotation",
            counts.mismatched
        );
        writeln!(err, "{line}")?;
        writeln!(report, "{line}")?;
    }
    if config.nmasking {
        let line = format!(
            "{} positions on chromosome {chr} were changed to 'N'",
            counts.new_n
        );
        writeln!(err, "{line}")?;
        writeln!(report, "{line}")?;
    }
    if config.full_sequence {
        let line = format!(
            "{} reference positions on chromosome {chr} were changed to the SNP alternative base\n",
            counts.new_snp
        );
        writeln!(err, "{line}")?;
        writeln!(report, "{line}")?;
    }
    writeln!(err)?;

    Ok((counts.new_n, counts.new_snp))
}

fn write_summary(
    config: &cli::Config,
    pass: Pass,
    label: &str,
    new_n_total: usize,
    new_snp_total: usize,
    report: &mut impl Write,
    err: &mut impl Write,
) -> Result<()> {
    let strain = &config.strain;
    let strain2 = config.strain2.as_deref().unwrap_or("");

    // Only the first pass leads with two blank lines, and only the third one drops the
    // separators from the report copy.
    let (lead, warn_who, report_who) = match pass {
        Pass::First => (
            "\n\n",
            format!("strain {strain}"),
            format!("strain {strain}"),
        ),
        Pass::Second => (
            "\n",
            format!("strain 2 [{label}]"),
            format!("strain 2 [{label}]"),
        ),
        Pass::DualHybrid => (
            "\n",
            format!("strain/strain 2 [{strain}/{strain2}]"),
            format!("strainstrain 2 [{strain}{strain2}]"),
        ),
    };

    if config.nmasking {
        writeln!(
            err,
            "{lead}Summary\n{new_n_total} Ns were newly introduced into the N-masked genome for {warn_who} in total"
        )?;
        writeln!(
            report,
            "\nSummary\n{new_n_total} Ns were newly introduced into the N-masked genome for {report_who} in total"
        )?;
    }
    if config.full_sequence {
        // The dual pass loses separators in both copies, but not the same ones: stderr says
        // "strainstrain 2 [A/B]" and the report says the same, while the N-masked line above
        // drops the slash inside the brackets too. Two adjacent lines, two different
        // omissions, reproduced exactly and raised upstream separately.
        let who = if pass == Pass::DualHybrid {
            format!("strainstrain 2 [{strain}/{strain2}]")
        } else {
            warn_who.clone()
        };
        let report_who = if pass == Pass::DualHybrid {
            format!("strainstrain 2 [{strain}/{strain2}]")
        } else {
            report_who.clone()
        };
        writeln!(
            err,
            "{new_snp_total} SNPs were newly introduced into the full sequence genome version for {who} in total\n"
        )?;
        writeln!(
            report,
            "{new_snp_total} SNPs were newly introduced into the full sequence genome version for {report_who} in total"
        )?;
    }
    // Only the first pass ends with a blank line: the Perl's `warn "\n"` sits in the strain-1
    // block of the main flow, and the two later blocks do not repeat it.
    if pass == Pass::First {
        writeln!(err)?;
    }
    Ok(())
}

fn summarise(config: &cli::Config, err: &mut impl Write) -> Result<()> {
    writeln!(err, "Summarising SNPsplit Genome Preparation Parameters")?;
    writeln!(err, "{}", "=".repeat(50))?;
    if !config.skip_filtering {
        writeln!(
            err,
            "Processing SNPs from VCF file:\t\t{}",
            config.vcf_file.as_ref().expect("checked").display()
        )?;
        writeln!(err, "Reading/filtering VCF file:\t\tYes (default)")?;
    } else {
        writeln!(err, "Reading/filtering VCF file:\t\tNo (skipped by user)")?;
    }
    writeln!(err, "Reference genome:\t\t\t{}", config.genome_folder)?;
    writeln!(
        err,
        "N-masking:\t\t\t\t{}",
        if config.nmasking { "Yes" } else { "No" }
    )?;
    writeln!(
        err,
        "Full SNP genome:\t\t\t{}",
        if config.full_sequence { "Yes" } else { "No" }
    )?;
    writeln!(err, "SNP strain:\t\t\t\t{}", config.strain)?;
    if let Some(strain2) = &config.strain2 {
        writeln!(err, "SNP strain 2:\t\t\t\t{strain2}")?;
    }
    if config.dual_hybrid {
        writeln!(
            err,
            "Dual hybrid, new Ref/SNP:\t\t{}/{}",
            config.strain,
            config.strain2.as_deref().unwrap_or("")
        )?;
    }
    writeln!(err)?;
    Ok(())
}

/// A whole-genome naming mismatch produces an unmodified copy of the reference rather than an
/// error, so it is caught up front rather than per chromosome.
fn check_chromosome_names_agree(
    genome: &genome::Genome,
    chroms: &BTreeSet<String>,
    err: &mut impl Write,
) -> Result<()> {
    if genome.keys().any(|c| chroms.contains(c)) {
        return Ok(());
    }

    writeln!(
        err,
        "\nNone of the chromosome names in the VCF file were found in the reference genome.\nA rather common mistake might be that the VCF file was downloaded from Ensembl (who use chromosome names such as 1, 2, X, MT)\nbut the genome from UCSC (who use chromosome names such as chr1, chr2, chrX, chrM)"
    )?;
    writeln!(err, "The chromosome names in the SNP annotation were:")?;
    for c in chroms {
        writeln!(err, "{c}")?;
    }
    writeln!(
        err,
        "The chromosome names in the reference genome folder were:"
    )?;
    for c in genome.keys() {
        writeln!(err, "{c}")?;
    }
    anyhow::bail!(
        "[FATAL ERROR] Please ensure that the same version of the genome is used for both VCF annotations and reference genome (FastA files). Exiting...\n\n"
    );
}
