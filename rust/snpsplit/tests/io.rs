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

use snpsplit::io::{sort_by_name, strnum_cmp};

/// The point of an external sort is the path where it spills, so this forces spilling with a
/// tiny in-memory budget rather than trusting the in-memory path to represent it.
#[test]
fn sorting_spills_to_disk_and_still_produces_one_ordered_stream() {
    let dir = TempDir::new().unwrap();
    let src = dir.path().join("unsorted.sam");
    let dst = dir.path().join("sorted.bam");

    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n");
    // Descending, so nothing is accidentally already in order.
    for i in (1..=200).rev() {
        sam.push_str(&format!(
            "read{i}\t0\tchr1\t{i}\t42\t4M\t*\t0\t0\tACGT\tIIII\n"
        ));
    }
    std::fs::File::create(&src)
        .unwrap()
        .write_all(sam.as_bytes())
        .unwrap();

    let header = RecordReader::open(&src).unwrap().header().clone();
    sort_by_name(&src, &dst, &header, 16).unwrap();

    let names = names_of(&dst);
    assert_eq!(names.len(), 200);

    let mut expected: Vec<String> = (1..=200).map(|i| format!("read{i}")).collect();
    expected.sort_by(|a, b| strnum_cmp(a.as_bytes(), b.as_bytes()));
    assert_eq!(names, expected);
    assert_eq!(names[0], "read1");
    assert_eq!(
        names[1],
        "read2",
        "numeric ordering lost: got {:?}",
        &names[..3]
    );
}

/// The comparator is checked against a handful of hand-picked names elsewhere. This checks
/// the whole sort, at a size that spills, against real samtools output. Skipped when
/// samtools is absent.
#[test]
fn the_whole_sort_matches_samtools_over_five_thousand_shuffled_names() {
    let dir = TempDir::new().unwrap();
    let src = dir.path().join("big.sam");
    let dst = dir.path().join("big.bam");

    // A fixed shuffle, so a failure is reproducible. Values are unrelated to their order.
    let mut order: Vec<usize> = (1..=5000).collect();
    let mut state: u64 = 0x2545F4914F6CDD1D;
    for i in (1..order.len()).rev() {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let j = (state >> 33) as usize % (i + 1);
        order.swap(i, j);
    }

    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n");
    for i in &order {
        sam.push_str(&format!(
            "read{i}\t0\tchr1\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\n"
        ));
    }
    std::fs::File::create(&src)
        .unwrap()
        .write_all(sam.as_bytes())
        .unwrap();

    let Ok(out) = std::process::Command::new("samtools")
        .args(["sort", "-n", "-O", "sam", src.to_str().unwrap()])
        .output()
    else {
        return;
    };
    if !out.status.success() {
        return;
    }
    let samtools_order: Vec<String> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(samtools_order.len(), 5000);

    let header = RecordReader::open(&src).unwrap().header().clone();
    sort_by_name(&src, &dst, &header, 500).unwrap();

    let ours = names_of(&dst);
    let first_divergence = ours
        .iter()
        .zip(&samtools_order)
        .position(|(a, b)| a != b)
        .map(|i| format!("at {i}: ours {} vs samtools {}", ours[i], samtools_order[i]));
    assert_eq!(
        first_divergence, None,
        "name sort diverged from samtools: {first_divergence:?}",
    );
    assert_eq!(ours, samtools_order);
}
