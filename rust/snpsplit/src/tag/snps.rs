//! The SNP annotation table.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};

/// Reference and alternative base at one position.
#[derive(Debug, Clone, Copy)]
pub struct Snp {
    pub reference: u8,
    pub alternative: u8,
}

/// Every SNP, by chromosome and 1-based position.
pub type Table = HashMap<String, HashMap<usize, Snp>>;

/// Read the SNP annotation, returning the table and how many positions it holds.
///
/// The format is the five-column track the genome preparation writes: id, chromosome,
/// position, strand, `Ref/SNP`. Lines that do not look like that are counted and skipped
/// rather than stored under an empty chromosome name, which is what used to happen.
pub fn read_snps(path: &Path, err: &mut impl Write) -> Result<(Table, usize)> {
    let file =
        File::open(path).with_context(|| format!("Failed to read from {}", path.display()))?;
    let reader: Box<dyn Read> = if path.to_string_lossy().ends_with(".gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    writeln!(
        err,
        "Storing SNP positions provided in '{}'",
        path.display()
    )?;

    let mut table: Table = HashMap::new();
    let mut count = 0usize;
    let mut skipped = 0usize;

    for line in BufReader::new(reader).lines() {
        let line = line?;
        let line: String = line.chars().filter(|c| *c != '\r' && *c != '\n').collect();
        if line.is_empty() {
            continue;
        }
        // The per-chromosome files written by the genome preparation open with a >chr header.
        if line.starts_with('>') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let Some(diff) = fields.get(4) else {
            skipped += 1;
            continue;
        };
        if fields.len() < 5 || !diff.contains('/') {
            skipped += 1;
            continue;
        }

        let (Some(chr), Some(pos)) = (fields.get(1), fields.get(2)) else {
            skipped += 1;
            continue;
        };
        let Ok(pos) = pos.parse::<usize>() else {
            skipped += 1;
            continue;
        };

        let (reference, alternative) = diff.split_once('/').expect("checked above");
        count += 1;
        table.entry(chr.to_string()).or_default().insert(
            pos,
            Snp {
                reference: *reference.as_bytes().first().unwrap_or(&b'N'),
                alternative: *alternative.as_bytes().first().unwrap_or(&b'N'),
            },
        );
    }

    writeln!(err, "Stored {count} positions in total\n")?;
    if skipped > 0 {
        writeln!(
            err,
            "Skipped {skipped} line(s) in '{}' that did not look like SNP annotations (expected: ID, chromosome, position, strand, Ref/SNP)\n",
            path.display()
        )?;
    }

    Ok((table, count))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn written(contents: &str) -> (TempDir, std::path::PathBuf) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("snps.txt");
        File::create(&path)
            .unwrap()
            .write_all(contents.as_bytes())
            .unwrap();
        (dir, path)
    }

    #[test]
    fn reads_the_five_column_track() {
        let (_dir, path) = written("1\tchr1\t30\t1\tT/A\n2\tchr1\t60\t1\tG/A\n");
        let (table, count) = read_snps(&path, &mut Vec::new()).unwrap();
        assert_eq!(count, 2);
        assert_eq!(table["chr1"][&30].reference, b'T');
        assert_eq!(table["chr1"][&30].alternative, b'A');
    }

    /// The per-chromosome files carry a `>name` header, which is not a SNP and not an error.
    #[test]
    fn a_chromosome_header_is_neither_stored_nor_counted_as_skipped() {
        let (_dir, path) = written(">1\n1\tchr1\t30\t1\tT/A\n");
        let mut log = Vec::new();
        let (_table, count) = read_snps(&path, &mut log).unwrap();
        assert_eq!(count, 1);
        assert!(!String::from_utf8_lossy(&log).contains("Skipped"));
    }

    /// A malformed line used to be stored under an empty chromosome name and counted as a
    /// SNP, which inflated the total and matched nothing.
    #[test]
    fn a_malformed_line_is_skipped_and_reported() {
        let (_dir, path) = written("1\tchr1\t30\t1\tT/A\nrubbish\n");
        let mut log = Vec::new();
        let (table, count) = read_snps(&path, &mut log).unwrap();
        assert_eq!(count, 1);
        assert_eq!(table.len(), 1);
        assert!(String::from_utf8_lossy(&log).contains("Skipped 1 line"));
    }
}
