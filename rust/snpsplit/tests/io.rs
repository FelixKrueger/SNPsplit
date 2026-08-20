//! The BAM/SAM layer, exercised against files shaped like the ones the fixtures carry.

use std::io::Write;
use std::path::PathBuf;

use snpsplit::io::RecordReader;
use tempfile::TempDir;

const SAM: &str = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:1000
read1\t0\tchr1\t1\t42\t4M\t*\t0\t0\tACGT\tIIII
read2\t16\tchr1\t10\t42\t4M\t*\t0\t0\tTTTT\tIIII
";

fn sam_file(dir: &TempDir) -> PathBuf {
    let path = dir.path().join("in.sam");
    std::fs::File::create(&path)
        .unwrap()
        .write_all(SAM.as_bytes())
        .unwrap();
    path
}

fn names_of(path: &std::path::Path) -> Vec<String> {
    RecordReader::open(path)
        .unwrap()
        .map(|r| String::from_utf8(r.unwrap().name().unwrap().to_vec()).unwrap())
        .collect()
}

#[test]
fn reads_the_header_and_every_record_from_sam() {
    let dir = TempDir::new().unwrap();
    let path = sam_file(&dir);
    let reader = RecordReader::open(&path).unwrap();

    assert_eq!(reader.header().reference_sequences().len(), 1);
    assert_eq!(names_of(&path), vec!["read1", "read2"]);
}

#[test]
fn a_missing_file_reports_the_path_the_perl_way() {
    let Err(err) = RecordReader::open(std::path::Path::new("/nonexistent/reads.bam")) else {
        panic!("opening a file that is not there should fail");
    };
    assert!(
        err.to_string().contains("/nonexistent/reads.bam"),
        "error did not name the file: {err}",
    );
}

/// A real fixture input, converted to BAM by samtools, read back through the noodles path.
/// Skipped rather than failed when samtools is absent: the fixture runners already require
/// it, so this is a convenience check, not the gate.
#[test]
fn reads_a_real_fixture_converted_to_bam() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../test/fixtures/se_basic/input.sam");
    if !fixture.exists() {
        return;
    }

    let dir = TempDir::new().unwrap();
    let bam = dir.path().join("input.bam");
    let converted = std::process::Command::new("samtools")
        .args(["view", "-bS", fixture.to_str().unwrap()])
        .stdout(std::fs::File::create(&bam).unwrap())
        .status();
    let Ok(status) = converted else { return };
    if !status.success() {
        return;
    }

    let from_sam = names_of(&fixture);
    let from_bam = names_of(&bam);
    assert!(!from_sam.is_empty(), "fixture had no records");
    assert_eq!(from_sam, from_bam, "BAM and SAM paths disagreed");
}
