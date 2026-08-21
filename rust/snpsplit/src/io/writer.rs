//! Writing alignments, to either SAM or BAM.

use std::fs::File;
use std::io::BufWriter;
use std::num::NonZero;
use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
// The record-writing method lives on this trait rather than on the concrete writers, which
// is what lets one enum cover both formats without a per-format record type.
use noodles_sam::alignment::io::Write as AlignmentWrite;
use noodles_sam::header::record::value::Map;
use noodles_sam::header::record::value::map::Program;
use noodles_sam::header::record::value::map::program::tag;

use super::Format;

/// Record this tool in the header's `@PG` chain.
///
/// Every program that writes an alignment file is supposed to say so, and the chain is how a
/// user reconstructs what touched their data. The Perl pipeline records only the samtools
/// invocations it shells out to, so the tool that actually did the work never appears; a
/// self-contained build has no excuse for that, and records itself.
///
/// `add` assigns a unique ID from the prefix and links `PP:` to the previous entry, which is
/// the same bookkeeping samtools does.
pub fn add_pg_line(
    header: &mut Header,
    id_prefix: &str,
    version: &str,
    command_line: &str,
) -> Result<()> {
    let program = Map::<Program>::builder()
        .insert(tag::NAME, id_prefix)
        .insert(tag::VERSION, version)
        .insert(tag::COMMAND_LINE, command_line)
        .build()
        .context("Failed to build the @PG record")?;

    header
        .programs_mut()
        .add(id_prefix, program)
        .context("Failed to add the @PG record to the header")?;
    Ok(())
}

/// Render one record as the SAM line it would be written as.
///
/// `--verbose` echoes every alignment, and the Perl can simply print the line it read because
/// it reads text. This build reads records, so the line has to be rendered back.
pub fn render_sam_line(header: &Header, record: &RecordBuf) -> Result<String> {
    let mut buffer: Vec<u8> = Vec::new();
    {
        let mut writer = noodles_sam::io::Writer::new(&mut buffer);
        writer.write_alignment_record(header, record)?;
    }
    Ok(String::from_utf8_lossy(&buffer).trim_end().to_string())
}

/// The command line as invoked, for the `CL:` field.
pub fn command_line() -> String {
    std::env::args().collect::<Vec<_>>().join(" ")
}

enum Inner {
    Sam(noodles_sam::io::Writer<BufWriter<File>>),
    Bam(noodles_bam::io::Writer<noodles_bgzf::io::MultithreadedWriter<File>>),
}

/// Writes alignment records, and reports failure rather than swallowing it.
///
/// `finish` is not optional and not a `Drop` impl: BGZF has a terminating block and both
/// forms have a final flush, and a `Drop` that cannot return an error is exactly how
/// tag2sort came to report success on a failed write (#116).
pub struct RecordWriter {
    inner: Inner,
}

impl RecordWriter {
    /// Create `path` for single-threaded writing.
    pub fn create(path: &Path, format: Format, header: &Header) -> Result<Self> {
        Self::create_with_workers(
            path,
            format,
            header,
            NonZero::new(1).expect("1 is not zero"),
        )
    }

    /// Create `path`, write `header` into it, and compress BGZF across `workers` threads.
    ///
    /// BGZF is a sequence of independently compressed blocks, so spreading the compression
    /// over threads cannot reorder anything: the blocks are emitted in the order their
    /// contents were written. The compressed bytes may differ from a single-threaded
    /// encoding; the decompressed content does not, and what the fixtures compare is
    /// `samtools view` output rather than BGZF bytes.
    pub fn create_with_workers(
        path: &Path,
        format: Format,
        header: &Header,
        workers: NonZero<usize>,
    ) -> Result<Self> {
        let creating = || format!("Failed to write to file '{}'", path.display());
        let file = File::create(path).with_context(creating)?;

        let inner = match format {
            Format::Sam => {
                let mut w = noodles_sam::io::Writer::new(BufWriter::new(file));
                w.write_header(header).with_context(creating)?;
                Inner::Sam(w)
            }
            Format::Bam => {
                let encoder =
                    noodles_bgzf::io::MultithreadedWriter::with_worker_count(workers, file);
                let mut w = noodles_bam::io::Writer::from(encoder);
                w.write_header(header).with_context(creating)?;
                Inner::Bam(w)
            }
        };

        Ok(Self { inner })
    }

    /// Append one record.
    pub fn write(&mut self, header: &Header, record: &RecordBuf) -> Result<()> {
        match &mut self.inner {
            Inner::Sam(w) => w.write_alignment_record(header, record)?,
            Inner::Bam(w) => w.write_alignment_record(header, record)?,
        }
        Ok(())
    }

    /// Flush and close, surfacing any error. Always call this.
    pub fn finish(self) -> Result<()> {
        match self.inner {
            Inner::Sam(w) => {
                let mut inner = w.into_inner();
                std::io::Write::flush(&mut inner).context("Failed to flush SAM output")?;
            }
            Inner::Bam(w) => {
                // The multithreaded encoder finishes on Drop, but Drop cannot report a
                // failure, which is the bug this whole type exists to avoid (#116).
                let mut encoder = w.into_inner();
                encoder.finish().context("Failed to finalise BAM output")?;
            }
        }
        Ok(())
    }
}
