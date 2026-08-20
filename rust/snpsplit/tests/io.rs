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

use snpsplit::io::{Format, RecordWriter};

#[test]
fn a_bam_round_trip_preserves_every_record() {
    let dir = TempDir::new().unwrap();
    let src = sam_file(&dir);
    let dst = dir.path().join("out.bam");

    let mut reader = RecordReader::open(&src).unwrap();
    let header = reader.header().clone();
    let mut writer = RecordWriter::create(&dst, Format::Bam, &header).unwrap();
    for record in reader.by_ref() {
        writer.write(&header, &record.unwrap()).unwrap();
    }
    writer.finish().unwrap();

    assert_eq!(names_of(&dst), vec!["read1", "read2"]);
}

/// tag2sort held BAM writers open for a whole run and did not check them (#116). A write
/// that cannot land has to surface rather than leaving a zero-byte file and a success
/// report.
#[test]
fn creating_a_bam_where_it_cannot_be_written_is_an_error_not_a_silent_success() {
    let dir = TempDir::new().unwrap();
    let src = sam_file(&dir);
    let readonly = dir.path().join("readonly");
    std::fs::create_dir(&readonly).unwrap();
    let mut perms = std::fs::metadata(&readonly).unwrap().permissions();
    perms.set_readonly(true);
    std::fs::set_permissions(&readonly, perms).unwrap();

    let header = RecordReader::open(&src).unwrap().header().clone();
    let result = RecordWriter::create(&readonly.join("out.bam"), Format::Bam, &header);

    assert!(
        result.is_err(),
        "creating a BAM in a read-only directory should fail",
    );
}

/// What the fixtures compare is samtools-rendered SAM text, so the contract for our BAM is
/// that samtools can read it. Skipped when samtools is absent.
#[test]
fn samtools_can_read_the_bam_we_write() {
    let dir = TempDir::new().unwrap();
    let src = sam_file(&dir);
    let dst = dir.path().join("out.bam");

    let mut reader = RecordReader::open(&src).unwrap();
    let header = reader.header().clone();
    let mut writer = RecordWriter::create(&dst, Format::Bam, &header).unwrap();
    for record in reader.by_ref() {
        writer.write(&header, &record.unwrap()).unwrap();
    }
    writer.finish().unwrap();

    let Ok(out) = std::process::Command::new("samtools")
        .args(["view", "-h", dst.to_str().unwrap()])
        .output()
    else {
        return;
    };
    if !out.status.success() {
        panic!(
            "samtools refused our BAM: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("read1"), "samtools output lost read1: {text}");
    assert!(text.contains("read2"), "samtools output lost read2: {text}");
}
