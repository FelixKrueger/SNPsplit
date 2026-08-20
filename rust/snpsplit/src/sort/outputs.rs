//! The output files a sorting run produces, and the names it derives for them.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

use crate::io::{Format, RecordWriter};

use super::cli::Options;

/// One destination: the records that land in it, and the name it was announced under.
pub struct Sink {
    // Taken by `finish`, which consumes the writer so its final flush can report a failure.
    writer: Option<RecordWriter>,
    header: Header,
}

impl Sink {
    fn create(path: &str, format: Format, header: &Header) -> Result<Self> {
        let writer = RecordWriter::create(Path::new(path), format, header)
            .with_context(|| format!("Unable to write to file '{path}'"))?;
        Ok(Self {
            writer: Some(writer),
            header: header.clone(),
        })
    }

    fn write(&mut self, record: &RecordBuf) -> Result<()> {
        match self.writer.as_mut() {
            Some(writer) => writer.write(&self.header, record),
            None => Ok(()),
        }
    }

    /// Close the writer, surfacing the failure a dropped one would swallow.
    fn finish(&mut self) -> Result<()> {
        match self.writer.take() {
            Some(writer) => writer.finish(),
            None => Ok(()),
        }
    }
}

/// Every file a run writes.
pub struct Outputs {
    pub genome1: Sink,
    pub genome2: Sink,
    pub unassigned: Sink,
    pub conflicting: Option<Sink>,
    pub genome1_st: Option<Sink>,
    pub genome2_st: Option<Sink>,
    pub unassigned_st: Option<Sink>,
    pub conflicting_st: Option<Sink>,
    pub g1_ua: Option<Sink>,
    pub g2_ua: Option<Sink>,
    pub g1_g2: Option<Sink>,
    pub report: BufWriter<File>,
    pub yaml: BufWriter<File>,
    yaml_lines: Vec<String>,
}

/// Derive every output name from the input, the way `open_output_filehandles` does.
///
/// The input is normalised to a `.bam` stem first, so a `.sam` input produces the same names
/// as a `.bam` one. That normalisation is what makes `--sam` produce sorted output at all.
fn stem(infile: &str) -> String {
    for suffix in [".allele_flagged.bam", ".allele_flagged.sam"] {
        if let Some(base) = infile.strip_suffix(suffix) {
            return format!("{base}.bam");
        }
    }
    infile.to_string()
}

fn named(stem: &str, suffix: &str, bam: bool) -> String {
    let base = stem.strip_suffix("bam").unwrap_or(stem);
    let extension = if bam { "bam" } else { "sam" };
    format!("{base}{suffix}.{extension}")
}

impl Outputs {
    pub fn open(
        opts: &Options,
        infile: &str,
        output_dir: &str,
        header: &Header,
        err: &mut impl Write,
    ) -> Result<Self> {
        let bam = !opts.sam;
        let format = if bam { Format::Bam } else { Format::Sam };
        let stem = stem(infile);

        let name = |suffix: &str| named(&stem, suffix, bam);
        let path = |n: &str| format!("{output_dir}{n}");

        let (genome1_file, genome2_file, unassigned_file) = if opts.hic {
            (name("G1_G1"), name("G2_G2"), name("UA_UA"))
        } else {
            (name("genome1"), name("genome2"), name("unassigned"))
        };
        let conflicting_file = name("conflicting");

        let base = stem.strip_suffix("bam").unwrap_or(&stem);
        let report_file = format!("{base}SNPsplit_sort.txt");
        let yaml_file = format!("{base}SNPsplit_sort.yaml");

        let mut outputs = Self {
            genome1: Sink::create(&path(&genome1_file), format, header)?,
            genome2: Sink::create(&path(&genome2_file), format, header)?,
            unassigned: Sink::create(&path(&unassigned_file), format, header)?,
            conflicting: opts
                .conflict
                .then(|| Sink::create(&path(&conflicting_file), format, header))
                .transpose()?,
            genome1_st: None,
            genome2_st: None,
            unassigned_st: None,
            conflicting_st: None,
            g1_ua: None,
            g2_ua: None,
            g1_g2: None,
            report: BufWriter::new(
                File::create(path(&report_file))
                    .with_context(|| format!("Unable to write to file '{report_file}'"))?,
            ),
            yaml: BufWriter::new(
                File::create(path(&yaml_file))
                    .with_context(|| format!("Unable to write to file '{yaml_file}'"))?,
            ),
            yaml_lines: Vec::new(),
        };

        let (g1_ua_file, g2_ua_file, g1_g2_file) = (name("G1_UA"), name("G2_UA"), name("G1_G2"));
        if opts.hic {
            outputs.g1_ua = Some(Sink::create(&path(&g1_ua_file), format, header)?);
            outputs.g2_ua = Some(Sink::create(&path(&g2_ua_file), format, header)?);
            outputs.g1_g2 = Some(Sink::create(&path(&g1_g2_file), format, header)?);
        }

        let (genome1_st, genome2_st, unassigned_st, conflicting_st) = (
            name("genome1_st"),
            name("genome2_st"),
            name("unassigned_st"),
            name("conflicting_st"),
        );
        if opts.singletons {
            outputs.genome1_st = Some(Sink::create(&path(&genome1_st), format, header)?);
            outputs.genome2_st = Some(Sink::create(&path(&genome2_st), format, header)?);
            outputs.unassigned_st = Some(Sink::create(&path(&unassigned_st), format, header)?);
            outputs.conflicting_st = opts
                .conflict
                .then(|| Sink::create(&path(&conflicting_st), format, header))
                .transpose()?;
        }

        writeln!(err, "Writing YAML report to {yaml_file}\n")?;

        // Collected first rather than written through a closure: the report and stderr both
        // need every line, and the one line that goes to stderr alone sits in the middle.
        let mut both: Vec<String> = vec![format!("Input file:\t\t\t\t\t\t'{infile}'")];
        // Announced to the terminal only: the report is the file being described, so naming
        // itself inside itself would be noise. It sits between the first line and the rest.
        let terminal_only = format!("Writing SNPsplit-sort report to:\t\t\t'{report_file}'");

        both.push(format!(
            "Writing unassigned reads to:\t\t\t\t'{unassigned_file}'"
        ));
        both.push(format!(
            "Writing genome 1-specific reads to:\t\t\t'{genome1_file}'"
        ));
        both.push(format!(
            "Writing genome 2-specific reads to:\t\t\t'{genome2_file}'"
        ));

        if opts.hic {
            both.push(format!("Writing G1/UA reads to:\t\t\t\t\t'{g1_ua_file}'"));
            both.push(format!("Writing G2/UA reads to:\t\t\t\t\t'{g2_ua_file}'"));
            both.push(format!("Writing G1/G2 reads to:\t\t\t\t\t'{g1_g2_file}'"));
        }

        if opts.conflict {
            both.push(format!(
                "Writing reads with conflicting number of SNPs to:\t'{conflicting_file}'\n"
            ));
        }

        if opts.singletons {
            both.push(format!(
                "Writing unassigned singleton reads to:\t\t\t'{unassigned_st}'"
            ));
            both.push(format!(
                "Writing genome 1-specific singleton reads to:\t\t'{genome1_st}'"
            ));
            both.push(format!(
                "Writing genome 2-specific singleton reads to:\t\t'{genome2_st}'"
            ));
            if opts.conflict {
                both.push(format!(
                    "Writing singleton reads with conflicting number of SNPs to:'{conflicting_st}'\n"
                ));
            }
        }

        for (index, line) in both.iter().enumerate() {
            writeln!(err, "{line}")?;
            writeln!(outputs.report, "{line}")?;
            if index == 0 {
                writeln!(err, "{terminal_only}")?;
            }
        }

        Ok(outputs)
    }

    pub fn write(sink: &mut Sink, record: &RecordBuf) -> Result<()> {
        sink.write(record)
    }

    /// Record one YAML entry, in the order they were produced.
    pub fn yaml(&mut self, key: &str, value: impl std::fmt::Display) {
        self.yaml_lines.push(format!("{key}: {value}"));
    }

    /// Close every record writer, surfacing any failure.
    pub fn finish(&mut self) -> Result<()> {
        let sinks = [
            Some(&mut self.genome1),
            Some(&mut self.genome2),
            Some(&mut self.unassigned),
            self.conflicting.as_mut(),
            self.genome1_st.as_mut(),
            self.genome2_st.as_mut(),
            self.unassigned_st.as_mut(),
            self.conflicting_st.as_mut(),
            self.g1_ua.as_mut(),
            self.g2_ua.as_mut(),
            self.g1_g2.as_mut(),
        ];
        for sink in sinks.into_iter().flatten() {
            sink.finish()?;
        }
        Ok(())
    }

    /// Flush the report and the YAML file.
    pub fn close_report(&mut self) -> Result<()> {
        for line in &self.yaml_lines {
            writeln!(self.yaml, "{line}")?;
        }
        self.yaml.flush()?;
        self.report.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_sam_input_produces_the_same_names_as_a_bam_one() {
        assert_eq!(stem("s.allele_flagged.bam"), "s.bam");
        assert_eq!(stem("s.allele_flagged.sam"), "s.bam");
    }

    #[test]
    fn an_input_that_is_not_allele_flagged_keeps_its_name_as_the_stem() {
        assert_eq!(stem("s.bam"), "s.bam");
    }

    #[test]
    fn output_names_follow_the_stem_and_the_output_format() {
        assert_eq!(named("s.bam", "genome1", true), "s.genome1.bam");
        assert_eq!(named("s.bam", "genome1", false), "s.genome1.sam");
        assert_eq!(named("s.bam", "G1_UA", true), "s.G1_UA.bam");
    }
}
