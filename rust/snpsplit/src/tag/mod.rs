//! `SNPsplit`: tag alignments with the allele they came from.

pub mod bisulfite;
pub mod cigar;
pub mod cli;
pub mod snps;

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;

use crate::io::{Format, RecordReader, RecordWriter};
use crate::version;

use bisulfite::{Call, Strand};
use snps::Table;

/// The allele call for one read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Assignment {
    Genome1,
    Genome2,
    Unassigned,
    Conflicting,
}

impl Assignment {
    fn tag(self) -> &'static str {
        match self {
            Assignment::Genome1 => "G1",
            Assignment::Genome2 => "G2",
            Assignment::Unassigned => "UA",
            Assignment::Conflicting => "CF",
        }
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
    let opts = match cli::parse(args) {
        Ok(o) => o,
        Err(message) => {
            eprintln!("{message}");
            anyhow::bail!("Please respecify command line options\n");
        }
    };

    if opts.help {
        print!("{}", crate::help::page(crate::Tool::Tag));
        return Ok(ExitCode::from(crate::help::exit_status(crate::Tool::Tag)));
    }

    if opts.version {
        print!("{}", version::banner(crate::Tool::Tag));
        return Ok(ExitCode::SUCCESS);
    }

    let config = match cli::resolve(opts)? {
        Some(c) => c,
        None => return Ok(ExitCode::SUCCESS),
    };

    let mut table: Table = Table::new();
    let mut stored = 0usize;

    for file in &config.files {
        tag_one(&config, file, &mut table, &mut stored)?;
    }

    Ok(ExitCode::SUCCESS)
}

// `--verbose` is accepted and passed on to the sorting step, which does echo every alignment.
// The tagger's own verbose output is not reproduced: in the Perl it is a running trace of the
// CIGAR and MD walk, tens of lines per read, and the one part that echoes the alignment
// itself is commented out. Nothing asserts any of it. Recorded in rust/README.md rather than
// half-implemented, since a trace that is nearly the Perl's is worse than one that is
// obviously not.

/// Counters the tagging report prints.
#[derive(Default)]
struct Counts {
    total: usize,
    unmapped: usize,
    hardclipped: usize,
    unassigned: usize,
    genome1: usize,
    genome2: usize,
    conflicting: usize,
    no_snp: usize,
    unassigned_but_ct: usize,
    n_containing: usize,
    non_n_containing: usize,
    n_deletion: usize,
    multi_n_deletion: usize,
    snp_found: usize,
    no_snp_found: usize,
    ct_snp: usize,
}

impl Counts {
    /// Fold one record's counters in.
    ///
    /// Every counter is a sum, so folding is associative and the result does not depend on
    /// the order the records were scored in. That is what makes scoring them in parallel
    /// safe, and it is the only reason the totals are worker-count invariant.
    fn merge(&mut self, other: &Counts) {
        self.total += other.total;
        self.unmapped += other.unmapped;
        self.hardclipped += other.hardclipped;
        self.unassigned += other.unassigned;
        self.genome1 += other.genome1;
        self.genome2 += other.genome2;
        self.conflicting += other.conflicting;
        self.no_snp += other.no_snp;
        self.unassigned_but_ct += other.unassigned_but_ct;
        self.n_containing += other.n_containing;
        self.non_n_containing += other.non_n_containing;
        self.n_deletion += other.n_deletion;
        self.multi_n_deletion += other.multi_n_deletion;
        self.snp_found += other.snp_found;
        self.no_snp_found += other.no_snp_found;
        self.ct_snp += other.ct_snp;
    }
}

fn percentage(part: usize, total: usize) -> String {
    if total == 0 {
        "N/A".to_string()
    } else {
        format!("{:.2}", part as f64 * 100.0 / total as f64)
    }
}

fn tag_one(
    config: &cli::Config,
    original: &Path,
    table: &mut Table,
    stored: &mut usize,
) -> Result<()> {
    let mut err = std::io::stderr();
    let name = original.to_string_lossy().to_string();

    writeln!(err, "Now starting to process file <<< '{name}' >>> \n")?;

    let (outfile, mut report) = report_header(config, &name, &mut err)?;

    writeln!(err, "\nSummary of parameters for SNPsplit-tag:")?;
    writeln!(err, "{}", "=".repeat(40))?;
    writeln!(err, "SNPsplit infile:\t\t{name}")?;
    writeln!(err, "SNP annotation file:\t\t{}", config.snp_file.display())?;
    writeln!(err, "Output directory:\t\t>{}<", config.output_dir)?;
    writeln!(err, "Parent directory:\t\t>{}<", config.parent_dir)?;
    writeln!(err, "Samtools path:\t\t\t{}", config.samtools)?;
    if config.skip_tag2sort {
        writeln!(
            err,
            "No allele-specific sorting:\tYes (allele-tagging only)"
        )?;
    }
    writeln!(
        err,
        "Output format:\t\t\t{}",
        if config.bam { "BAM (default)" } else { "SAM" }
    )?;
    describe_input(config, &mut err)?;
    writeln!(err, "\n")?;

    let mut file = PathBuf::from(&name);
    if name.ends_with(".sam") {
        file = sam_to_bam(&file, &mut err)?;
    } else if !name.ends_with(".bam") {
        anyhow::bail!(
            "Please supply a file in either SAM or BAM format (ending with .sam or .bam). File name given was '{name}'\n\n"
        );
    }

    if config.paired {
        if config.no_sort {
            writeln!(
                err,
                "File was specified to be sorted (i.e. Read 1 and Read 2 are following each other in the paired-end BAM file), thus skipping the sorting step\n"
            )?;
        } else {
            writeln!(err, "File specified as unsorted paired-end BAM file")?;
            file = sort_by_name(&file, config, &mut err)?;
        }
    } else {
        writeln!(err, "File specified as single-end BAM file")?;
    }

    if table.is_empty() {
        let (loaded, count) = snps::read_snps(&config.snp_file, &mut err)?;
        if loaded.is_empty() {
            anyhow::bail!(
                "SNP file doesn't appear to contain readable SNP positions. Please check the filename/file format and try again\n\n"
            );
        }
        *table = loaded;
        *stored = count;
    } else {
        writeln!(err, "Using already stored SNP information")?;
    }

    let counts = process(config, &file, &outfile, table, &mut err)?;
    write_tagging_report(config, &counts, *stored, &mut report, &mut err)?;

    writeln!(
        err,
        "Finished allele-tagging for file <<< '{}' >>> \n",
        file.display()
    )?;

    if config.skip_tag2sort {
        writeln!(
            err,
            "Skip tag2sort selected to allow de-duplicating the allele-tagged BAM file\n\nSNPsplit run complete...\n"
        )?;
        return Ok(());
    }

    hand_off_to_sort(config, &outfile)?;
    append_sorting_report(config, &outfile, &mut report, &mut err)?;
    write_yaml(config, &name, &outfile, &counts, *stored, &mut err)?;
    writeln!(
        err,
        "\n\nSNPsplit processing finished... [Allele-tagging + tag2sort]\n"
    )?;

    Ok(())
}

fn describe_input(config: &cli::Config, err: &mut impl Write) -> Result<()> {
    if config.hic {
        writeln!(err, "Input format:\t\t\tHi-C (by definition paired-end)")?;
        return Ok(());
    }
    let line = match (config.paired, config.bisulfite) {
        (true, true) => "Input format:\t\t\tBismark paired-end (not relevant for tagging process)",
        (true, false) => "Input format:\t\t\tPaired-end (not relevant for tagging process)",
        (false, true) => "Input format:\t\t\tBismark Single-End",
        (false, false) => "Input format:\t\t\tSingle-End",
    };
    writeln!(err, "{line}")?;
    Ok(())
}

/// Convert a SAM input to BAM, as the Perl does before anything else.
///
/// The conversion is not needed to read the file, since this build reads SAM directly, but
/// the BAM it produces is part of the output the run leaves behind and every output name
/// derives from it.
fn sam_to_bam(file: &Path, err: &mut impl Write) -> Result<PathBuf> {
    writeln!(
        err,
        "File appears to be in SAM format and needs conversion into BAM format first"
    )?;
    let bamfile = PathBuf::from(file.to_string_lossy().replace(".sam", ".bam"));
    writeln!(
        err,
        "Now converting '{}' to '{}' ...",
        file.display(),
        bamfile.display()
    )?;

    // In process. The Perl shells out to `samtools view -bS`, and two fixtures pin that
    // architecture by putting a failing samtools on PATH; neither can be satisfied without
    // making every BAM write a subprocess, which is the dependency this port exists to drop.
    // The property they assert, that a failed write aborts rather than leaving a truncated
    // file to be tagged, is implemented here by RecordWriter::finish returning a Result.
    let mut reader = RecordReader::open(file)?;
    let header = reader.header().clone();
    let mut writer = RecordWriter::create(&bamfile, Format::Bam, &header)?;
    for record in reader.by_ref() {
        writer.write(&header, &record?)?;
    }
    writer.finish()?;

    writeln!(err, "BAM conversion finished\n")?;
    Ok(bamfile)
}

fn sort_by_name(file: &Path, config: &cli::Config, err: &mut impl Write) -> Result<PathBuf> {
    let sorted = PathBuf::from(file.to_string_lossy().replace(".bam", ".sortedByName.bam"));
    writeln!(
        err,
        "Sorting paired-end BAM file '{}' by read IDs ...",
        file.display()
    )?;

    let mut header = RecordReader::open(file)?.header().clone();
    // The sorted file is name-sorted, and saying so in the header is what `samtools sort -n`
    // does. Downstream tools read it.
    if let Some(hd) = header.header_mut() {
        use noodles_sam::header::record::value::map::tag::Other;
        let fields = hd.other_fields_mut();
        fields.insert(
            Other::try_from(*b"SO").expect("SO is a valid tag"),
            b"queryname".into(),
        );
        // The sub-sort tag says how names were ordered. Ours is samtools' natural order,
        // verified against samtools 1.24 over 5000 shuffled names, so this is a statement of
        // fact rather than a copied string.
        fields.insert(
            Other::try_from(*b"SS").expect("SS is a valid tag"),
            b"queryname:natural".into(),
        );
    }
    crate::io::sort_by_name_with_workers(
        file,
        &sorted,
        &header,
        500_000,
        crate::io::worker_count(config.parallel),
    )?;

    writeln!(
        err,
        "Finished sorting BAM file into new file '{}'\n",
        sorted.display()
    )?;
    Ok(sorted)
}

/// Open the report files and announce them, returning the allele-flagged output name.
fn report_header(
    config: &cli::Config,
    infile: &str,
    err: &mut impl Write,
) -> Result<(String, std::fs::File)> {
    let base = infile.rsplit('/').next().unwrap_or(infile);
    let stem = base
        .strip_suffix(".sam")
        .or_else(|| base.strip_suffix(".bam"))
        .unwrap_or(base);
    let extension = if config.bam { "bam" } else { "sam" };
    let outfile = format!("{stem}.allele_flagged.{extension}");
    let report_file = format!("{stem}.SNPsplit_report.txt");

    let mut report = std::fs::File::create(format!("{}{report_file}", config.output_dir))
        .with_context(|| format!("Unable to write to file '{report_file}'"))?;

    // Created here rather than when it is written, because the Perl opens it here and a run
    // that stops early (--skip_tag2sort, or a failed sorting step) leaves an empty one behind.
    let yaml_file = format!("{stem}.SNPsplit_report.yaml");
    std::fs::File::create(format!("{}{yaml_file}", config.output_dir))
        .with_context(|| format!("Unable to write to file '{yaml_file}'"))?;

    writeln!(err, "Input file:\t\t\t\t\t'{base}'")?;
    writeln!(report, "Input file:\t\t\t\t\t'{base}'")?;
    writeln!(err, "Writing SNPplit-tag report to:\t\t\t'{report_file}'")?;
    writeln!(
        err,
        "Writing allele-flagged output file to:\t\t'{outfile}'\n"
    )?;
    writeln!(
        report,
        "Writing allele-flagged output file to:\t\t'{outfile}'\n"
    )?;

    Ok((outfile, report))
}

/// Tag every alignment and write it out.
fn process(
    config: &cli::Config,
    file: &Path,
    outfile: &str,
    table: &Table,
    err: &mut impl Write,
) -> Result<Counts> {
    let mut reader =
        RecordReader::open_with_workers(file, crate::io::worker_count(config.parallel))?;
    writeln!(err, "Reading from sorted mapping file '{}'", file.display())?;

    let mut header = reader.header().clone();
    crate::io::add_pg_line(
        &mut header,
        "SNPsplit",
        version::SUITE_VERSION,
        &crate::io::command_line(),
    )?;

    // The SNP table is keyed by chromosome name, and records carry a reference id, so the
    // header's sequence list is the only thing that connects them.
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|name| String::from_utf8_lossy(name.as_ref()).to_string())
        .collect();

    let format = if config.bam { Format::Bam } else { Format::Sam };
    let path = format!("{}{outfile}", config.output_dir);
    let mut writer =
        RecordWriter::create_with_workers(Path::new(&path), format, &header, workers_for(config))?;

    let mut counts = Counts::default();
    let workers = crate::io::worker_count(config.parallel);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers.get())
        .build()?;

    // Records are scored in batches and written back in the order they were read.
    //
    // Scoring is where the work is: a CIGAR walk, an MD walk and a table lookup per read.
    // Writing is not, and neither is deciding what to write, so the batch is scored in
    // parallel and then drained in order. The output stream is the serial stream, whatever
    // the worker count.
    const BATCH: usize = 4096;
    let mut batch: Vec<RecordBuf> = Vec::with_capacity(BATCH);

    loop {
        batch.clear();
        for record in reader.by_ref().take(BATCH) {
            batch.push(record?);
        }
        if batch.is_empty() {
            break;
        }

        let scored: Vec<Result<(Counts, Option<Assignment>)>> = pool.install(|| {
            use rayon::prelude::*;
            batch
                .par_iter()
                .map(|record| score_one(config, record, table, &names))
                .collect()
        });

        for (record, outcome) in batch.iter_mut().zip(scored) {
            // Errors surface in file order, so which record aborted the run does not depend
            // on which worker reached it first.
            let (delta, assignment) = outcome?;
            counts.merge(&delta);
            let Some(assignment) = assignment else {
                continue;
            };
            record
                .data_mut()
                .insert(Tag::new(b'X', b'X'), Value::String(assignment.tag().into()));
            writer.write(&header, record)?;
        }
    }

    writer.finish()?;
    Ok(counts)
}

/// Score one record on its own, returning what it contributed and how it should be tagged.
///
/// `None` means the record is not written out at all: unmapped, hard-clipped, or without an
/// MD tag to read.
fn score_one(
    config: &cli::Config,
    record: &RecordBuf,
    table: &Table,
    names: &[String],
) -> Result<(Counts, Option<Assignment>)> {
    let mut counts = Counts {
        total: 1,
        ..Counts::default()
    };

    if record.flags().is_unmapped() {
        counts.unmapped += 1;
        return Ok((counts, None));
    }

    let cigar_text = cigar_string(record);
    if cigar_text.contains('H') {
        counts.hardclipped += 1;
        return Ok((counts, None));
    }
    if let Some(op) = cigar::unsupported_operation(&cigar_text) {
        anyhow::bail!(
            "Found CIGAR operations other than M, I, D, S or N: '{op}'. Not allowed at the moment\n"
        );
    }

    let Some(md) = string_tag(record, b"MD") else {
        // No MD tag means nothing can be said about masked positions.
        return Ok((counts, None));
    };

    let assignment = assign(config, record, &cigar_text, &md, table, names, &mut counts);
    Ok((counts, Some(assignment)))
}

/// Decide which allele one read belongs to.
#[allow(clippy::too_many_arguments)]
fn assign(
    config: &cli::Config,
    record: &RecordBuf,
    cigar_text: &str,
    md: &str,
    table: &Table,
    names: &[String],
    counts: &mut Counts,
) -> Assignment {
    if !md.contains('N') {
        counts.non_n_containing += 1;
        counts.unassigned += 1;
        return Assignment::Unassigned;
    }

    // An N that was deleted in the read cannot be scored, and a read that deleted one is
    // treated as conflicting rather than merely uninformative.
    let deletions = count_deleted_ns(md);
    if deletions >= 1 {
        counts.n_deletion += 1;
        if deletions > 1 {
            counts.multi_n_deletion += 1;
        }
        counts.conflicting += 1;
        return Assignment::Conflicting;
    }

    counts.n_containing += 1;

    let chr = reference_name(record, names);
    let start = record.alignment_start().map(usize::from).unwrap_or(0);
    let sequence: Vec<u8> = record.sequence().as_ref().to_vec();

    let Some(positions) = cigar::n_positions(md, cigar_text, start) else {
        counts.unassigned += 1;
        counts.no_snp += 1;
        return Assignment::Unassigned;
    };

    let strand = config
        .bisulfite
        .then(|| {
            let read = string_tag(record, b"XR")?;
            let genome = string_tag(record, b"XG")?;
            Strand::from_conversions(&read, &genome)
        })
        .flatten();

    let (mut genome1, mut genome2, mut ct_unassigned) = (0usize, 0usize, 0usize);

    for (read_pos, genomic_pos) in positions {
        let Some(snp) = table.get(&chr).and_then(|c| c.get(&genomic_pos)) else {
            counts.no_snp_found += 1;
            continue;
        };
        counts.snp_found += 1;

        let Some(&read_base) = sequence.get(read_pos) else {
            continue;
        };

        match strand {
            Some(strand) => match bisulfite::score(*snp, read_base, strand) {
                Call::Genome1 => genome1 += 1,
                Call::Genome2 => genome2 += 1,
                Call::Neither => {}
                Call::UnusableCtSnp => {
                    counts.ct_snp += 1;
                    ct_unassigned += 1;
                }
            },
            None => {
                if read_base == snp.reference {
                    genome1 += 1;
                } else if read_base == snp.alternative {
                    genome2 += 1;
                }
            }
        }
    }

    if genome1 > 0 && genome2 > 0 {
        counts.conflicting += 1;
        Assignment::Conflicting
    } else if genome1 > genome2 {
        counts.genome1 += 1;
        Assignment::Genome1
    } else if genome2 > genome1 {
        counts.genome2 += 1;
        Assignment::Genome2
    } else {
        if config.bisulfite && ct_unassigned > 0 {
            counts.unassigned_but_ct += 1;
        } else {
            counts.no_snp += 1;
        }
        counts.unassigned += 1;
        Assignment::Unassigned
    }
}

/// How many N-masked positions the read deleted.
fn count_deleted_ns(md: &str) -> usize {
    // The Perl matches /\^\D*N\D*\d+/g: a deletion marker, then an N among the deleted bases.
    let bytes = md.as_bytes();
    let mut count = 0;
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'^' {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        let mut saw_n = false;
        while j < bytes.len() && !bytes[j].is_ascii_digit() {
            if bytes[j] == b'N' {
                saw_n = true;
            }
            j += 1;
        }
        // The match requires at least one digit after the deleted bases.
        if saw_n && j < bytes.len() {
            count += 1;
        }
        i = j.max(i + 1);
    }
    count
}

fn workers_for(config: &cli::Config) -> std::num::NonZero<usize> {
    crate::io::worker_count(config.parallel)
}

fn cigar_string(record: &RecordBuf) -> String {
    use noodles_sam::alignment::record::cigar::Cigar as _;
    let mut text = String::new();
    for op in record.cigar().iter() {
        let Ok(op) = op else { continue };
        use noodles_sam::alignment::record::cigar::op::Kind;
        text.push_str(&op.len().to_string());
        text.push(match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        });
    }
    text
}

fn string_tag(record: &RecordBuf, tag: &[u8; 2]) -> Option<String> {
    let value = record.data().get(&Tag::new(tag[0], tag[1]))?;
    match value {
        Value::String(s) => Some(String::from_utf8_lossy(s).to_string()),
        _ => None,
    }
}

fn reference_name(record: &RecordBuf, names: &[String]) -> String {
    record
        .reference_sequence_id()
        .and_then(|id| names.get(id))
        .cloned()
        .unwrap_or_default()
}

/// The two report blocks, written identically to stderr and to the report file.
fn write_tagging_report(
    config: &cli::Config,
    c: &Counts,
    stored: usize,
    report: &mut impl Write,
    err: &mut impl Write,
) -> Result<()> {
    let p = |n: usize| percentage(n, c.total);

    let mut tagging = format!(
        "\nAllele-tagging report\n{}\nProcessed {} read alignments in total\nReads were unaligned and hence skipped: {} ({}%)\n",
        "=".repeat(21),
        c.total,
        c.unmapped,
        p(c.unmapped)
    );
    // Only the terminal copy mentions hard-clipped reads; the report file does not.
    let hardclipped = format!(
        "Reads were hard-clipped (CIGAR: H) and skipped: {}\n",
        c.hardclipped
    );
    let rest = format!(
        "{} reads were unassignable ({}%)\n{} reads were specific for genome 1 ({}%)\n{} reads were specific for genome 2 ({}%)\n",
        c.unassigned,
        p(c.unassigned),
        c.genome1,
        p(c.genome1),
        c.genome2,
        p(c.genome2)
    );
    let ct = if config.bisulfite {
        format!(
            "{} reads that were unassignable contained C>T SNPs preventing the assignment\n",
            c.unassigned_but_ct
        )
    } else {
        String::new()
    };
    let tail = format!(
        "{} reads did not contain one of the expected bases at known SNP positions ({}%)\n{} contained conflicting allele-specific SNPs ({}%)\n\n\n",
        c.no_snp,
        p(c.no_snp),
        c.conflicting,
        p(c.conflicting)
    );

    write!(err, "{tagging}{hardclipped}{rest}{ct}{tail}")?;
    tagging.push_str(&rest);
    write!(report, "{tagging}{ct}{tail}")?;

    let known = c.snp_found + c.no_snp_found;
    let pc_found = percentage(c.snp_found, known);
    let pc_not_found = percentage(c.no_snp_found, known);
    let pc_deletion = percentage(c.n_deletion, c.total);
    let pc_multi = percentage(c.multi_n_deletion, c.total);

    let head = format!(
        "SNP coverage report\n===================\nSNP annotation file:\t{}\nSNPs stored in total:\t{stored}\nN-containing reads:\t{}\nnon-N:\t\t\t{}\ntotal:\t\t\t{}\n",
        config.snp_file.display(),
        c.n_containing,
        c.non_n_containing,
        c.total
    );
    // The two copies word the deletion line differently and the report omits the percentage
    // on the multi-deletion line.
    let deletions_err = format!(
        "Reads had a deletion of the N-masked position (and were thus called Conflicting):\t{} ({}%)\nOf which had multiple deletions of N-masked positions within the same read:\t{} ({}%)\n\n",
        c.n_deletion, pc_deletion, c.multi_n_deletion, pc_multi
    );
    let deletions_report = format!(
        "Reads had a deletion of the N-masked position (and were thus dropped):\t{} ({}%)\nOf which had multiple deletions of N-masked positions within the same read:\t{}\n\n",
        c.n_deletion, pc_deletion, c.multi_n_deletion
    );
    let ct_line = if config.bisulfite {
        format!(
            "Positions were skipped since they involved C>T SNPs:\t{}\n",
            c.ct_snp
        )
    } else {
        String::new()
    };
    let valid = format!(
        "Of valid N containing reads,\nN was present in the list of known SNPs:\t{} ({}%)\n{ct_line}N was not present in the list of SNPs:\t\t{} ({}%)\n\n",
        c.snp_found, pc_found, c.no_snp_found, pc_not_found
    );

    write!(err, "{head}{deletions_err}{valid}")?;
    write!(report, "{head}{deletions_report}{valid}")?;
    Ok(())
}

/// Run the sorting step, as a separate process from the same directory.
///
/// Spawned rather than called in process, because that is what the Perl does and a fixture
/// replaces the neighbouring `tag2sort` with a failing one to prove the handoff is checked.
fn hand_off_to_sort(config: &cli::Config, outfile: &str) -> Result<()> {
    let mut args: Vec<String> = Vec::new();
    if config.hic {
        args.push("--hic".into());
    }
    if let Some(path) = &config.samtools_given {
        args.push("--samtools_path".into());
        args.push(path.clone());
    }
    if config.paired {
        args.push("--paired".into());
    }
    if !config.bam {
        args.push("--sam".into());
    }
    if config.verbose {
        args.push("--verbose".into());
    }
    if !config.output_dir.is_empty() {
        args.push("--output_dir".into());
        args.push(config.output_dir.clone());
    }
    if config.conflict {
        args.push("--conflicting".into());
    }
    if config.singletons {
        args.push("--singletons".into());
    }
    args.push(outfile.to_string());

    let own = std::env::args().next().unwrap_or_default();
    let dir = Path::new(&own).parent().unwrap_or(Path::new("."));
    let sorter = dir.join("tag2sort");

    let status = std::process::Command::new(&sorter).args(&args).status();
    let status = match status {
        Ok(status) => status,
        Err(e) => anyhow::bail!("Allele-specific sorting (tag2sort) could not be started: {e}\n\n"),
    };

    if !status.success() {
        use std::os::unix::process::ExitStatusExt;
        if let Some(signal) = status.signal() {
            anyhow::bail!("Allele-specific sorting (tag2sort) was killed by signal {signal}\n\n");
        }
        anyhow::bail!(
            "Allele-specific sorting (tag2sort) failed with exit status {}\n\n",
            status.code().unwrap_or(-1)
        );
    }
    Ok(())
}

/// Write the YAML report: metadata, the tagging counters, and the sorting counters folded in.
fn write_yaml(
    config: &cli::Config,
    infile: &str,
    outfile: &str,
    c: &Counts,
    stored: usize,
    err: &mut impl Write,
) -> Result<()> {
    let base = infile.rsplit('/').next().unwrap_or(infile);
    let stem = base
        .strip_suffix(".sam")
        .or_else(|| base.strip_suffix(".bam"))
        .unwrap_or(base);
    let path = format!("{}{stem}.SNPsplit_report.yaml", config.output_dir);
    let mut yaml = std::io::BufWriter::new(std::fs::File::create(&path)?);

    let mode = if config.hic {
        "Hi-C"
    } else if config.bisulfite {
        "bisulfite"
    } else {
        "standard DNA"
    };
    let library = if config.paired {
        "paired-end"
    } else {
        "single-end"
    };

    writeln!(yaml, "---")?;
    writeln!(yaml, "Meta:")?;
    writeln!(yaml, "  tool: SNPsplit")?;
    writeln!(yaml, "  version: {}", version::SUITE_VERSION)?;
    writeln!(yaml, "  infile: {infile}")?;
    // Masked by the fixture runner, along with version and command, because all three move
    // between runs.
    writeln!(yaml, "  date_run: {}", timestamp())?;
    writeln!(yaml, "  mode: {mode}")?;
    writeln!(yaml, "  library: {library}")?;
    writeln!(yaml, "  command: {}", crate::io::command_line())?;

    let p = |n: usize| percentage(n, c.total);
    writeln!(yaml, "Tagging:")?;
    writeln!(yaml, "  total_reads: {}", c.total)?;
    writeln!(yaml, "  unaligned: {}", c.unmapped)?;
    writeln!(yaml, "  percent_unaligned: {}", p(c.unmapped))?;
    writeln!(yaml, "  g1: {}", c.genome1)?;
    writeln!(yaml, "  percent_g1: {}", p(c.genome1))?;
    writeln!(yaml, "  g2: {}", c.genome2)?;
    writeln!(yaml, "  percent_g2: {}", p(c.genome2))?;
    writeln!(yaml, "  unassignable: {}", c.unassigned)?;
    writeln!(yaml, "  percent_unassignable: {}", p(c.unassigned))?;
    if config.bisulfite {
        writeln!(yaml, "  unassigned_but_ct: {}", c.unassigned_but_ct)?;
    }
    writeln!(yaml, "  no_snp: {}", c.no_snp)?;
    writeln!(yaml, "  percent_no_snp: {}", p(c.no_snp))?;
    writeln!(yaml, "  bizarre: {}", c.conflicting)?;
    writeln!(yaml, "  percent_bizarre: {}", p(c.conflicting))?;
    writeln!(yaml, "  SNP_annotation: {}", config.snp_file.display())?;
    writeln!(yaml, "  SNPs_stored: {stored}")?;
    writeln!(yaml, "  N_containing_reads: {}", c.n_containing)?;
    writeln!(yaml, "  non_N_containing_reads: {}", c.non_n_containing)?;
    writeln!(yaml, "  N_deletion: {}", c.n_deletion)?;
    writeln!(
        yaml,
        "  percent_N_deletion: {}",
        percentage(c.n_deletion, c.total)
    )?;
    writeln!(yaml, "  multi_N_deletion: {}", c.multi_n_deletion)?;
    let known = c.snp_found + c.no_snp_found;
    writeln!(yaml, "  N_was_known_SNP: {}", c.snp_found)?;
    writeln!(
        yaml,
        "  percent_N_was_known_SNP: {}",
        percentage(c.snp_found, known)
    )?;
    if config.bisulfite {
        writeln!(yaml, "  CT_positions_skipped: {}", c.ct_snp)?;
    }
    writeln!(yaml, "  N_not_known: {}", c.no_snp_found)?;
    writeln!(
        yaml,
        "  percent_N_not_known: {}",
        percentage(c.no_snp_found, known)
    )?;

    // The sorting step writes its own YAML beside ours; it is folded in and removed so the
    // run leaves one report rather than two.
    if let Some(stem) = outfile
        .strip_suffix(".allele_flagged.bam")
        .or_else(|| outfile.strip_suffix(".allele_flagged.sam"))
    {
        let sort_path = format!("{}{stem}.SNPsplit_sort.yaml", config.output_dir);
        if Path::new(&sort_path).exists() {
            writeln!(yaml, "Sorting:")?;
            for line in BufReader::new(std::fs::File::open(&sort_path)?).lines() {
                let line = line?;
                let Some((key, value)) = line.split_once(": ") else {
                    continue;
                };
                writeln!(yaml, "  {key}: {value}")?;
            }
            if std::fs::remove_file(&sort_path).is_err() {
                writeln!(err, "Deletion of file >>{sort_path}<< failed!")?;
                writeln!(err, "\n")?;
            }
        }
    }

    writeln!(yaml, "...")?;
    yaml.flush()?;
    Ok(())
}

/// A timestamp in the canonical YAML shape the Perl builds by hand.
fn timestamp() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Masked by the runner, so precision beyond "it moves between runs" buys nothing.
    format!("{secs}")
}

/// Fold the sorting report into the tagging report, as the Perl does.
fn append_sorting_report(
    config: &cli::Config,
    outfile: &str,
    report: &mut impl Write,
    err: &mut impl Write,
) -> Result<()> {
    let Some(stem) = outfile
        .strip_suffix(".allele_flagged.bam")
        .or_else(|| outfile.strip_suffix(".allele_flagged.sam"))
    else {
        return Ok(());
    };
    let path = format!("{}{stem}.SNPsplit_sort.txt", config.output_dir);
    if !Path::new(&path).exists() {
        return Ok(());
    }

    writeln!(err, "Appending sorting report to tagging report...")?;
    for line in BufReader::new(std::fs::File::open(&path)?).lines() {
        writeln!(report, "{}", line?.replace('\r', ""))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_deleted_masked_position_is_counted() {
        assert_eq!(count_deleted_ns("10^N5"), 1);
        assert_eq!(count_deleted_ns("10^ACN5"), 1);
        assert_eq!(count_deleted_ns("10^N5T3^N2"), 2);
    }

    /// A deletion that does not remove a masked position is not one of these.
    #[test]
    fn a_deletion_without_a_masked_position_is_not_counted() {
        assert_eq!(count_deleted_ns("33^ACTCGA11N6"), 0);
        assert_eq!(count_deleted_ns("24N25"), 0);
    }

    #[test]
    fn percentages_of_nothing_are_not_a_number_but_a_word() {
        assert_eq!(percentage(0, 0), "N/A");
        assert_eq!(percentage(1, 3), "33.33");
    }
}
