//! Writing alignments, to either SAM or BAM.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
// The record-writing method lives on this trait rather than on the concrete writers, which
// is what lets one enum cover both formats without a per-format record type.
use noodles_sam::alignment::io::Write as AlignmentWrite;

use super::Format;

enum Inner {
    Sam(noodles_sam::io::Writer<BufWriter<File>>),
    Bam(noodles_bam::io::Writer<noodles_bgzf::io::Writer<File>>),
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
    /// Create `path` and write `header` into it.
    pub fn create(path: &Path, format: Format, header: &Header) -> Result<Self> {
        let creating = || format!("Failed to write to file '{}'", path.display());
        let file = File::create(path).with_context(creating)?;

        let inner = match format {
            Format::Sam => {
                let mut w = noodles_sam::io::Writer::new(BufWriter::new(file));
                w.write_header(header).with_context(creating)?;
                Inner::Sam(w)
            }
            Format::Bam => {
                let mut w = noodles_bam::io::Writer::new(file);
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
            Inner::Bam(mut w) => {
                w.try_finish().context("Failed to finalise BAM output")?;
            }
        }
        Ok(())
    }
}
