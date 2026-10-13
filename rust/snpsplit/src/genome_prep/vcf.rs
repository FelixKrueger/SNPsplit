//! Reading and filtering the VCF: chromosome detection, strain detection, and the
//! homozygous high-confidence SNP filter that produces the `SNPs_<strain>/` track.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};

/// Open a VCF, transparently decompressing when the name ends in `.gz`.
///
/// The Perl forks `gunzip -c`. Decompressing in-process is what removes the runtime
/// dependency, and it is also why a path containing a space needs no special handling here.
pub fn open_vcf(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to read Input VCF file '{}'", path.display()))?;

    let name = path.to_string_lossy();
    let reader: Box<dyn Read> = if name.ends_with("gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

/// Strip the line terminators the Perl removes with `s/(\r|\n)//g`.
fn strip_eol(line: &str) -> &str {
    line.trim_end_matches(['\r', '\n'])
}

/// The chromosomes named by `##contig` header lines.
///
/// Reading stops at the first non-header line, as the Perl does: the contig declarations are
/// all in the header, and a whole-file scan of a mouse VCF would be minutes of work for
/// nothing.
pub fn detect_chroms(path: &Path) -> Result<Vec<String>> {
    let mut seen: Vec<String> = Vec::new();
    let reader = open_vcf(path)?;

    for line in reader.lines() {
        let line = line?;
        let line = strip_eol(&line);
        if !line.starts_with('#') {
            break;
        }
        if let Some(rest) = line.strip_prefix("##contig")
            && let Some(id) = extract_contig_id(rest)
            && !seen.iter().any(|c| c == &id)
        {
            seen.push(id);
        }
    }

    Ok(seen)
}

/// The `ID=...` value of a `##contig` line, up to the next comma.
///
/// Non-greedy, matching the Perl: a contig line may carry several commas and the first one
/// ends the identifier.
fn extract_contig_id(rest: &str) -> Option<String> {
    let start = rest.find("ID=")? + 3;
    let tail = &rest[start..];
    let end = tail.find(',')?;
    Some(tail[..end].to_string())
}

/// The strain names in the `#CHROM` header, mapped to their column index.
///
/// The first nine columns are fixed VCF fields; everything after them is a sample.
pub fn detect_strains(path: &Path) -> Result<BTreeMap<String, usize>> {
    let mut strains = BTreeMap::new();
    let reader = open_vcf(path)?;

    for line in reader.lines() {
        let line = line?;
        let line = strip_eol(&line);
        if line.starts_with("##") {
            continue;
        }
        if !line.starts_with('#') {
            break;
        }
        for (index, field) in line.split('\t').enumerate() {
            if index <= 8 {
                continue;
            }
            strains.insert(field.to_string(), index);
        }
    }

    Ok(strains)
}

/// Counters the filter reports, in the order the report prints them.
#[derive(Default)]
pub struct FilterReport {
    pub total: usize,
    pub indels: usize,
    pub homozygous: usize,
    pub high_confidence: usize,
    pub same_as_reference: usize,
    pub no_clear_alternative: usize,
    pub other_genotype: usize,
    pub low_confidence: usize,
}

/// A position that was homozygous in a strain but did not pass the filter, kept so the dual
/// hybrid comparison can tell "absent" from "present but low confidence".
#[derive(Default, Clone)]
pub struct Blacklisted {
    pub strain1_filter: Option<String>,
    pub strain2_filter: Option<String>,
}

/// Which of the two strains a filtering pass is for.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum StrainIdentity {
    First,
    Second,
}

/// Every high-confidence SNP, for the gzipped all-SNP list: by chromosome, then by position,
/// with the first call at a position winning, as the Perl's `%all_SNPs` does.
///
/// An MGP strain has around twenty million of them. A map entry per SNP, with the chromosome
/// name and the formatted line each in their own small allocation, came to about 160 bytes a
/// SNP, and the allocator does not hand that many small blocks back to the system, so it
/// stayed resident through the genome build that follows. One vector per chromosome, sorted
/// once at the end, is 32 bytes a SNP in a handful of large allocations.
#[derive(Default)]
struct AllSnps {
    by_chr: BTreeMap<String, Vec<SnpRecord>>,
}

struct SnpRecord {
    pos: u64,
    line: SnpLine,
}

enum SnpLine {
    /// The usual case: single-base alleles, and a position that prints as it parses.
    Compact { total: u64, reference: u8, alt: u8 },
    /// Anything else keeps its line verbatim, so the output never depends on the encoding.
    Verbatim(Box<str>),
}

impl AllSnps {
    fn insert(&mut self, chr: &str, pos: &str, total: usize, reference: &str, alt: &str) {
        let parsed = pos.parse::<u64>().unwrap_or(0);
        let canonical = !pos.is_empty()
            && pos.bytes().all(|b| b.is_ascii_digit())
            && (pos == "0" || !pos.starts_with('0'));
        let line = match (reference.as_bytes(), alt.as_bytes()) {
            ([r], [a]) if canonical => SnpLine::Compact {
                total: total as u64,
                reference: *r,
                alt: *a,
            },
            _ => SnpLine::Verbatim(
                format!("{total}\t{chr}\t{pos}\t1\t{reference}/{alt}").into_boxed_str(),
            ),
        };
        let record = SnpRecord { pos: parsed, line };
        match self.by_chr.get_mut(chr) {
            Some(records) => records.push(record),
            None => {
                self.by_chr.insert(chr.to_string(), vec![record]);
            }
        }
    }

    /// Order each chromosome by position and keep the first call at each. The sort is stable,
    /// so among calls at one position the earliest read comes first.
    fn finish(&mut self) {
        for records in self.by_chr.values_mut() {
            records.sort_by_key(|r| r.pos);
            records.dedup_by_key(|r| r.pos);
        }
    }

    fn write_lines(&self, out: &mut impl Write) -> std::io::Result<()> {
        for (chr, records) in &self.by_chr {
            for record in records {
                match &record.line {
                    SnpLine::Compact {
                        total,
                        reference,
                        alt,
                    } => writeln!(
                        out,
                        "{total}\t{chr}\t{}\t1\t{}/{}",
                        record.pos, *reference as char, *alt as char
                    )?,
                    SnpLine::Verbatim(line) => writeln!(out, "{line}")?,
                }
            }
        }
        Ok(())
    }
}

/// Everything the filter produces besides the files it writes.
pub struct FilterOutcome {
    pub report: FilterReport,
    pub all_snps_name: String,
}

/// Filter the VCF down to homozygous high-confidence SNPs for one strain.
///
/// Writes `SNPs_<strain>/chr<name>.txt` per chromosome, a `<strain>_SNP_filtering_report.txt`,
/// and a gzipped list of every SNP found.
#[allow(clippy::too_many_arguments)]
pub fn filter_snps(
    vcf: &Path,
    strain: &str,
    strain_index: usize,
    identity: StrainIdentity,
    chroms: &[String],
    genome_build: &str,
    v7: bool,
    dual_hybrid: bool,
    blacklist: &mut BTreeMap<String, Blacklisted>,
    out: &mut impl Write,
    err: &mut impl Write,
) -> Result<FilterOutcome> {
    let dir = format!("SNPs_{strain}");
    if !Path::new(&dir).is_dir() {
        writeln!(
            err,
            "Folder '{dir}' doesn't exist. Creating it for you...\n"
        )?;
        std::fs::create_dir(&dir)
            .with_context(|| format!("Failed to created directory {dir}\n"))?;
    }

    let mut handles: BTreeMap<String, std::io::BufWriter<File>> = BTreeMap::new();
    for chr in chroms {
        let filename = format!("SNPs_{strain}/chr{chr}.txt");
        let mut fh = std::io::BufWriter::new(File::create(&filename)?);
        writeln!(fh, ">{chr}")?;
        handles.insert(chr.clone(), fh);
    }

    let mut report = FilterReport::default();
    let mut all_snps = AllSnps::default();

    let mut format_index: Option<usize> = None;
    let mut info_index: Option<usize> = None;
    let mut gt_index: Option<usize> = None;
    let mut fi_index: Option<usize> = None;

    let reader = open_vcf(vcf)?;
    let mut input_line = 0usize;
    for line in reader.lines() {
        let line = line?;
        input_line += 1;
        let line = strip_eol(&line);

        if line.starts_with("##") {
            continue;
        }

        if line.starts_with("#CHROM") {
            let fields: Vec<&str> = line.split('\t').collect();
            let name = fields.get(strain_index).copied().unwrap_or("");
            writeln!(err, "Analysing SNP fields for name >{name}<")?;

            for (index, field) in fields.iter().enumerate() {
                if *field == "FORMAT" {
                    format_index = Some(index);
                }
                if *field == "INFO" {
                    info_index = Some(index);
                }
            }

            match format_index {
                Some(i) => writeln!(err, "Using FORMAT field index: {i}")?,
                None => {
                    anyhow::bail!(
                        "Failed to extract index of field 'FORMAT' from the VCF header (line {input_line}). Hmmm...\n"
                    );
                }
            }
            match info_index {
                Some(i) => writeln!(err, "Using INFO field index: {i}")?,
                None => {
                    anyhow::bail!(
                        "Failed to extract index of field 'INFO' from the VCF header (line {input_line}). Hmmm...\n"
                    );
                }
            }
            continue;
        }

        report.total += 1;
        if report.total.is_multiple_of(1_000_000) {
            writeln!(err, "processed {} lines", report.total)?;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let chr = get(0);
        let pos = get(1);
        let reference = get(3);
        let mut alt = get(4).to_string();
        let info = info_index.map(get).unwrap_or("");
        let format = format_index.map(get).unwrap_or("");
        let sample = get(strain_index);

        if gt_index.is_none() {
            writeln!(err, "GT index not defined, checking...")?;
            let keys: Vec<&str> = format.split(':').collect();
            match keys.iter().position(|k| *k == "GT") {
                Some(i) => {
                    gt_index = Some(i);
                    writeln!(err, "Setting GT index to >>{i}<<")?;
                }
                None => {
                    anyhow::bail!(
                        "Failed to extract GT index from the FORMAT field of the first SNP call (VCF line {input_line}). I am afraid this will need some manual looking into...\n"
                    );
                }
            }
            match keys.iter().position(|k| *k == "FI") {
                Some(i) => {
                    fi_index = Some(i);
                    writeln!(err, "Setting FI index to >>{i}<<")?;
                }
                None => {
                    anyhow::bail!(
                        "Failed to extract FI index from the FORMAT field of the first SNP call (VCF line {input_line}). I am afraid this will need some manual looking into...\n"
                    );
                }
            }
        }

        let sample_fields: Vec<&str> = sample.split(':').collect();
        let pick = |i: Option<usize>| i.and_then(|i| sample_fields.get(i).copied()).unwrap_or("");
        let gt = pick(gt_index);
        let fi = pick(fi_index);

        if v7 && info.contains("INDEL") {
            report.indels += 1;
            continue;
        }

        if reference
            .chars()
            .any(|c| !matches!(c, 'A' | 'T' | 'C' | 'G'))
        {
            writeln!(err, "ref was: {reference}; skipping")?;
            continue;
        }

        if alt
            .chars()
            .any(|c| !matches!(c, 'A' | 'T' | 'C' | 'G' | ','))
        {
            report.no_clear_alternative += 1;
            continue;
        }

        let location = format!("{chr}:{pos}");

        if dual_hybrid && fi != "1" {
            let entry = blacklist.entry(location.clone()).or_default();
            match identity {
                StrainIdentity::First => entry.strain1_filter = Some(fi.to_string()),
                StrainIdentity::Second => entry.strain2_filter = Some(fi.to_string()),
            }
        }

        // Only 1/1, 2/2 and 3/3 are homozygous alternative calls. The index into a
        // comma-separated ALT list is the genotype number minus one.
        let alt_index = match gt {
            "0/0" => {
                report.same_as_reference += 1;
                continue;
            }
            "1/1" => Some(0usize),
            "2/2" => Some(1usize),
            "3/3" => Some(2usize),
            _ => {
                report.other_genotype += 1;
                continue;
            }
        };

        report.homozygous += 1;
        if let Some(index) = alt_index
            && alt.contains(',')
        {
            alt = alt.split(',').nth(index).unwrap_or("").to_string();
        }

        if fi == "1" {
            report.high_confidence += 1;
            all_snps.insert(chr, pos, report.total, reference, &alt);
        } else {
            report.low_confidence += 1;
            continue;
        }

        // One handle exists per ##contig header line, so a chromosome the header never
        // declared has none. Checked here, at the point of use, so a low-confidence call on an
        // undeclared contig is still filtered out rather than aborting the run (#130).
        let Some(fh) = handles.get_mut(chr) else {
            anyhow::bail!(
                "The SNP call on VCF line {input_line} is on chromosome '{chr}', which is not declared in a ##contig line of the VCF header.\nPlease make sure the VCF header declares every chromosome its SNP calls use. Exiting...\n\n"
            );
        };
        writeln!(
            fh,
            "{}\t{chr}\t{pos}\t1\t{reference}/{alt}\t{sample}",
            report.total
        )?;
    }

    for (_, mut fh) in handles {
        fh.flush()?;
    }

    all_snps.finish();
    let all_snps_name = format!("all_SNPs_{strain}_{genome_build}.txt.gz");
    write_reports(
        strain,
        genome_build,
        &report,
        v7,
        &all_snps_name,
        &all_snps,
        out,
        err,
    )?;

    Ok(FilterOutcome {
        report,
        all_snps_name,
    })
}

#[allow(clippy::too_many_arguments)]
fn write_reports(
    strain: &str,
    genome_build: &str,
    report: &FilterReport,
    v7: bool,
    all_snps_name: &str,
    all_snps: &AllSnps,
    _out: &mut impl Write,
    err: &mut impl Write,
) -> Result<()> {
    let mut summary = String::new();
    summary.push_str(&format!(
        "SNP position summary for strain {strain} (based on genome build {genome_build})\n"
    ));
    summary.push_str(&"=".repeat(75));
    summary.push_str("\n\n");
    summary.push_str(&format!("Positions read in total:\t{}\n\n", report.total));
    if v7 {
        summary.push_str(&format!(
            "{}\tPositions were INDELs (and hence skipped)\n",
            report.indels
        ));
    }
    summary.push_str(&format!(
        "{}\tSNP were homozygous. Of these:\n",
        report.homozygous
    ));
    summary.push_str(&format!(
        "{}\tSNP were homozygous and passed high confidence filters and were thus included into the {strain} genome\n",
        report.high_confidence
    ));
    summary.push_str(&format!("\nNot included into {strain} genome:\n"));
    summary.push_str(&format!(
        "{}\thad the same sequence as the reference\n",
        report.same_as_reference
    ));
    summary.push_str(&format!(
        "{}\t\thad no clearly defined alternative base\n",
        report.no_clear_alternative
    ));
    summary.push_str(&format!(
        "{}\t\tCalls were neither 0/0 (same as reference) or 1/1, 2/2, 3/3 (homozygous SNP)\n",
        report.other_genotype
    ));
    summary.push_str(&format!(
        "{}\t\twere homozygous but the filtering call was low confidence\n\n",
        report.low_confidence
    ));

    write!(err, "\n{summary}")?;

    let report_name = format!("{strain}_SNP_filtering_report.txt");
    let mut fh = File::create(&report_name)
        .with_context(|| format!("Failed to write to file {report_name}"))?;
    write!(fh, "{summary}")?;

    writeln!(
        err,
        "Now printing a single list of all SNPs to >{all_snps_name}<..."
    )?;
    writeln!(
        fh,
        "Printed a single list of all SNPs to >{all_snps_name}<..."
    )?;
    drop(fh);

    if Path::new(all_snps_name).exists() {
        writeln!(
            err,
            "File '{all_snps_name}' existed in the folder already, overwriting it...\n"
        )?;
    }

    write_gzipped(all_snps_name, all_snps)?;

    writeln!(err, "complete\n")?;
    Ok(())
}

/// Write the all-SNP list, compressed.
///
/// Through `gzip -c` when there is one to run, in process when there is not.
///
/// Compressing in process would be the obvious choice, and the bytes it produces decompress
/// to exactly the same text. The reason to prefer the subprocess is that a *failing* gzip is
/// observable behaviour: SNPsplit reads this file, and a truncated one means it loads no SNPs
/// at all, so the Perl treats a gzip failure as fatal and a fixture pins that. With no
/// subprocess there is nothing to fail in that way, and the fixture could not express what it
/// tests.
///
/// The in-process fallback keeps the binary self-sufficient where the Perl would simply not
/// run. The only difference between the two paths is which failures are possible, never the
/// content.
fn write_gzipped(path: &str, snps: &AllSnps) -> Result<()> {
    use std::process::{Command, Stdio};

    let file = File::create(path).with_context(|| format!("Failed to write to file {path}"))?;

    let spawned = Command::new("gzip")
        .arg("-c")
        .stdin(Stdio::piped())
        .stdout(Stdio::from(file.try_clone()?))
        .spawn();

    let Ok(mut child) = spawned else {
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        snps.write_lines(&mut gz)?;
        gz.finish()
            .with_context(|| format!("Failed to write the SNP list to {path}"))?;
        return Ok(());
    };

    {
        // Buffered, since the lines are now formatted one write at a time. A gzip that has
        // already died closes the pipe, and the write fails rather than killing us with
        // SIGPIPE; either way the status check below is what reports it.
        let mut stdin = std::io::BufWriter::new(child.stdin.take().expect("stdin was piped"));
        let _ = snps.write_lines(&mut stdin).and_then(|()| stdin.flush());
    }

    let status = child.wait()?;
    if !status.success() {
        anyhow::bail!("Failed to write the SNP list to {path} via gzip: the file is incomplete\n");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_contig_id_stops_at_the_first_comma() {
        assert_eq!(
            extract_contig_id("=<ID=1,length=320>").as_deref(),
            Some("1")
        );
        assert_eq!(
            extract_contig_id("=<ID=chrX,length=1,assembly=x>").as_deref(),
            Some("chrX")
        );
    }

    /// The Perl's regex needs the comma; a contig line without one matches nothing. Reproduced
    /// rather than improved, because a contig that vanishes here changes which chromosomes are
    /// processed and that is observable.
    #[test]
    fn a_contig_id_without_a_trailing_comma_is_not_matched() {
        assert_eq!(extract_contig_id("=<ID=1>"), None);
    }

    #[test]
    fn eol_stripping_handles_windows_line_endings() {
        assert_eq!(strip_eol("a\tb\r\n"), "a\tb");
        assert_eq!(strip_eol("a\tb\n"), "a\tb");
        assert_eq!(strip_eol("a\tb"), "a\tb");
    }
}
