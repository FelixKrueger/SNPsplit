//! Reading alignments, from either SAM or BAM, as owned records.

use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result};
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;

/// Which of the two forms a file is in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Sam,
    Bam,
}

/// Decide the format from the file's first bytes, not its name.
///
/// The Perl decides by extension, which is why `--no_sorting` rejects anything not ending
/// `.bam`. Misnamed files are common enough downstream of other tools that reading the
/// magic is worth the two extra syscalls; the extension is not consulted at all.
pub fn sniff(path: &Path) -> Result<Format> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to open file '{}'", path.display()))?;
    let mut magic = [0u8; 2];
    let read = read_up_to(&mut file, &mut magic)?;
    file.seek(SeekFrom::Start(0))?;

    // gzip magic. BAM is BGZF, which is gzip; a plain-gzipped SAM is not something the Perl
    // accepts either, so treating gzip as BAM matches the reachable behaviour.
    if read == 2 && magic == [0x1f, 0x8b] {
        return Ok(Format::Bam);
    }
    Ok(Format::Sam)
}

fn read_up_to(file: &mut File, buf: &mut [u8]) -> Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match file.read(&mut buf[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok(filled)
}

enum Inner {
    Sam(noodles_sam::io::Reader<BufReader<File>>),
    Bam(noodles_bam::io::Reader<noodles_bgzf::io::Reader<File>>),
}

/// Streams owned alignment records out of a SAM or BAM file.
///
/// Records are owned (`RecordBuf`) rather than borrowed: the tagger holds a read while it
/// walks its CIGAR against the SNP table, and the sorter buffers mates, so a borrowed
/// record tied to the reader's buffer would not survive either.
pub struct RecordReader {
    inner: Inner,
    header: Header,
}

impl RecordReader {
    /// Open `path`, reading its header. Format is sniffed, not taken from the extension.
    pub fn open(path: &Path) -> Result<Self> {
        let opening = || format!("Failed to open file '{}'", path.display());

        let (inner, header) = match sniff(path)? {
            Format::Sam => {
                let file = File::open(path).with_context(opening)?;
                let mut reader = noodles_sam::io::Reader::new(BufReader::new(file));
                let header = reader.read_header().with_context(opening)?;
                (Inner::Sam(reader), header)
            }
            Format::Bam => {
                let file = File::open(path).with_context(opening)?;
                let mut reader = noodles_bam::io::Reader::new(file);
                let header = reader.read_header().with_context(opening)?;
                (Inner::Bam(reader), header)
            }
        };

        Ok(Self { inner, header })
    }

    /// The header as read from the file.
    pub fn header(&self) -> &Header {
        &self.header
    }
}

impl Iterator for RecordReader {
    type Item = Result<RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = RecordBuf::default();
        let read = match &mut self.inner {
            Inner::Sam(r) => r.read_record_buf(&self.header, &mut record),
            Inner::Bam(r) => r.read_record_buf(&self.header, &mut record),
        };
        match read {
            Ok(0) => None,
            Ok(_) => Some(Ok(record)),
            Err(e) => Some(Err(e.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn write(dir: &TempDir, name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = dir.path().join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(bytes)
            .unwrap();
        path
    }

    #[test]
    fn a_bgzf_magic_is_bam_whatever_the_name_says() {
        let dir = TempDir::new().unwrap();
        // BGZF is a gzip member with an extra BC field: 1f 8b 08 04 ...
        let p = write(&dir, "misnamed.sam", &[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0]);
        assert_eq!(sniff(&p).unwrap(), Format::Bam);
    }

    #[test]
    fn a_header_line_is_sam_whatever_the_name_says() {
        let dir = TempDir::new().unwrap();
        let p = write(&dir, "misnamed.bam", b"@HD\tVN:1.6\n");
        assert_eq!(sniff(&p).unwrap(), Format::Sam);
    }

    /// A headerless SAM is legal, so an empty or non-magic first two bytes must not be
    /// treated as BAM.
    #[test]
    fn a_headerless_record_line_is_sam() {
        let dir = TempDir::new().unwrap();
        let p = write(
            &dir,
            "reads.sam",
            b"read1\t0\tchr1\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\n",
        );
        assert_eq!(sniff(&p).unwrap(), Format::Sam);
    }

    #[test]
    fn an_empty_file_is_sam_so_the_reader_reports_no_records_rather_than_corruption() {
        let dir = TempDir::new().unwrap();
        let p = write(&dir, "empty.sam", b"");
        assert_eq!(sniff(&p).unwrap(), Format::Sam);
    }
}
