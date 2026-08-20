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

/// Reproduce the way Perl reports a `die` whose message has no trailing newline.
///
/// Perl appends ` at <script> line <n>, <IN> line <m>.` in that case, and four fixtures diff
/// it. The script line number is masked by the runner, so the value carries no information;
/// the Perl line is used anyway, as provenance. The input line number is *not* masked and is
/// genuinely useful, which is the half of this worth keeping.
///
/// The underlying defect is the missing newline: an internal message should not be leaking a
/// script line number to a user at all. Raised upstream separately; emulated here because
/// byte-identity is the gate and the fixture pins it.
fn perl_die(message: &str, perl_line: u32, input_line: usize) -> anyhow::Error {
    let script = std::env::args().next().unwrap_or_default();
    anyhow::anyhow!("{message} at {script} line {perl_line}, <IN> line {input_line}.\n")
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
    // Keyed by chr:pos, first occurrence wins, as the Perl's %all_SNPs does.
    let mut all_snps: BTreeMap<(String, u64), String> = BTreeMap::new();

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
                    return Err(perl_die(
                        "Failed to extract index of field 'FORMAT'. Hmmm...",
                        1001,
                        input_line,
                    ));
                }
            }
            match info_index {
                Some(i) => writeln!(err, "Using INFO field index: {i}")?,
                None => {
                    return Err(perl_die(
                        "Failed to extract index of field 'INFO'. Hmmm...",
                        1007,
                        input_line,
                    ));
                }
            }
            continue;
        }

        report.total += 1;
        if report.total % 1_000_000 == 0 {
            writeln!(err, "processed {} lines", report.total)?;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let chr = get(0);
        let pos = get(1);
        let reference = get(3);
        let mut alt = get(4).to_string();
        let info = info_index.map(&get).unwrap_or("");
        let format = format_index.map(&get).unwrap_or("");
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
                    return Err(perl_die(
                        "Failed to extract GT index. I am afraid this will need some manual looking into...",
                        1044,
                        input_line,
                    ));
                }
            }
            match keys.iter().position(|k| *k == "FI") {
                Some(i) => {
                    fi_index = Some(i);
                    writeln!(err, "Setting FI index to >>{i}<<")?;
                }
                None => {
                    return Err(perl_die(
                        "Failed to extract FI index. I am afraid this will need some manual looking into...",
                        1051,
                        input_line,
                    ));
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
            let snp = format!("{}\t{chr}\t{pos}\t1\t{reference}/{alt}", report.total);
            let key = (chr.to_string(), pos.parse::<u64>().unwrap_or(0));
            all_snps.entry(key).or_insert(snp);
        } else {
            report.low_confidence += 1;
            continue;
        }

        // A chromosome that appears in the body but was never declared in a ##contig header
        // has no filehandle. The Perl reaches `print {$fhs{$chr}}` with an undefined value and
        // the interpreter aborts; the fixture pins that abort, so the run has to end here too.
        let Some(fh) = handles.get_mut(chr) else {
            return Err(perl_die(
                "Can't use an undefined value as a symbol reference",
                1190,
                input_line,
            ));
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
    all_snps: &BTreeMap<(String, u64), String>,
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

    let file = File::create(all_snps_name)
        .with_context(|| format!("Failed to write to file {all_snps_name}"))?;
    let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    for snp in all_snps.values() {
        writeln!(gz, "{snp}")?;
    }
    gz.finish()
        .with_context(|| format!("Failed to write the SNP list to {all_snps_name}"))?;

    writeln!(err, "complete\n")?;
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
