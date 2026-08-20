//! Sorting alignments by read name, the way `samtools sort -n` does.
//!
//! The Perl delegates this to samtools (`SNPsplit:1394`). Doing it here is what lets the
//! binary drop its samtools runtime dependency, but it means reproducing samtools'
//! comparator exactly: output order is part of the byte-identity contract, and samtools
//! does not sort names byte-wise.

use std::cmp::Ordering;

/// samtools' `strnum_cmp`: alternating text and numeric segments.
///
/// This is a port of samtools' implementation, not a reimplementation of the idea. Two
/// behaviours here are easy to get wrong and both are observable:
///
/// - digit runs compare by value, so `read2` precedes `read10`;
/// - when two runs have the same value but different numbers of leading zeroes, the one
///   with more leading zeroes sorts first, so `read007` precedes `read7`. They are not
///   equal. Verified against `samtools sort -n` (1.24), not assumed.
pub fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (mut ia, mut ib) = (0usize, 0usize);

    while ia < a.len() && ib < b.len() {
        if a[ia].is_ascii_digit() && b[ib].is_ascii_digit() {
            // Leading zeroes do not contribute to the value, but how many were skipped is
            // remembered in ia/ib and breaks the tie further down.
            while ia < a.len() && a[ia] == b'0' {
                ia += 1;
            }
            while ib < b.len() && b[ib] == b'0' {
                ib += 1;
            }

            while ia < a.len()
                && ib < b.len()
                && a[ia].is_ascii_digit()
                && b[ib].is_ascii_digit()
                && a[ia] == b[ib]
            {
                ia += 1;
                ib += 1;
            }

            let a_digit = ia < a.len() && a[ia].is_ascii_digit();
            let b_digit = ib < b.len() && b[ib].is_ascii_digit();

            if a_digit && b_digit {
                // Both runs continue and differ here. Whichever run is longer is the larger
                // number, whatever the digits at this position say.
                let mut i = 0;
                while ia + i < a.len()
                    && ib + i < b.len()
                    && a[ia + i].is_ascii_digit()
                    && b[ib + i].is_ascii_digit()
                {
                    i += 1;
                }
                return if ia + i < a.len() && a[ia + i].is_ascii_digit() {
                    Ordering::Greater
                } else if ib + i < b.len() && b[ib + i].is_ascii_digit() {
                    Ordering::Less
                } else {
                    a[ia].cmp(&b[ib])
                };
            } else if a_digit {
                return Ordering::Greater;
            } else if b_digit {
                return Ordering::Less;
            } else if ia != ib {
                // Same value, different zero padding. More padding sorts first.
                return if ia < ib {
                    Ordering::Greater
                } else {
                    Ordering::Less
                };
            }
        } else {
            if a[ia] != b[ib] {
                return a[ia].cmp(&b[ib]);
            }
            ia += 1;
            ib += 1;
        }
    }

    // Whichever still has bytes left is the longer, hence greater, name.
    if ia < a.len() {
        Ordering::Greater
    } else if ib < b.len() {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

use std::collections::BinaryHeap;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::Result;
use noodles_sam::Header;
use noodles_sam::alignment::RecordBuf;
use rayon::slice::ParallelSliceMut;

use super::{Format, RecordReader, RecordWriter};

/// Sort `src` by read name into `dst`, single-threaded.
pub fn sort_by_name(
    src: &Path,
    dst: &Path,
    header: &Header,
    max_records_in_memory: usize,
) -> Result<()> {
    sort_by_name_with_workers(
        src,
        dst,
        header,
        max_records_in_memory,
        NonZero::new(1).expect("1 is not zero"),
    )
}

/// Sort `src` by read name into `dst`, using `workers` threads.
///
/// Runs of at most `max_records_in_memory` records are sorted in memory and spilled to
/// temporary BAM files, then merged. `samtools sort` uses the same shape. The budget is a
/// record count rather than a byte count because record size varies little within a run and
/// counting bytes would mean measuring every record twice.
///
/// Two things take the worker count: BGZF compression of the spill files and of the output,
/// and the in-memory sort of each run. The merge does not, and will not. The merge is where
/// the final order is decided, and a faster wrong order is worse than a slower right one.
///
/// The run sort is `par_sort_by`, which is stable, matching the `sort_by` it replaces. An
/// unstable parallel sort would reorder records that compare equal, and how many compared
/// equal within a run depends on where the run boundaries fell, so the output would start
/// depending on the memory budget. Worker-count invariance is not free; this is one of the
/// places it is paid for.
pub fn sort_by_name_with_workers(
    src: &Path,
    dst: &Path,
    header: &Header,
    max_records_in_memory: usize,
    workers: NonZero<usize>,
) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers.get())
        .build()?;
    if matches!(super::sniff(src)?, Format::Bam) {
        return sort_bam_raw(&pool, src, dst, header, max_records_in_memory, workers);
    }

    let spill_dir = tempfile::tempdir()?;
    let mut runs: Vec<PathBuf> = Vec::new();
    let mut buffer: Vec<RecordBuf> = Vec::with_capacity(max_records_in_memory);

    for record in RecordReader::open_with_workers(src, workers)? {
        buffer.push(record?);
        if buffer.len() >= max_records_in_memory {
            runs.push(spill(
                &pool,
                &mut buffer,
                header,
                spill_dir.path(),
                runs.len(),
                workers,
            )?);
        }
    }

    // Nothing spilled means the whole input fitted, which is the common case at fixture
    // size: sort and write it straight out rather than paying for a one-way merge.
    if runs.is_empty() {
        pool.install(|| buffer.par_sort_by(compare));
        let mut writer = RecordWriter::create_with_workers(dst, Format::Bam, header, workers)?;
        for record in &buffer {
            writer.write(header, record)?;
        }
        return writer.finish();
    }

    if !buffer.is_empty() {
        runs.push(spill(
            &pool,
            &mut buffer,
            header,
            spill_dir.path(),
            runs.len(),
            workers,
        )?);
    }

    merge(&runs, dst, header, workers)
}

fn spill(
    pool: &rayon::ThreadPool,
    buffer: &mut Vec<RecordBuf>,
    header: &Header,
    dir: &Path,
    index: usize,
    workers: NonZero<usize>,
) -> Result<PathBuf> {
    pool.install(|| buffer.par_sort_by(compare));
    let path = dir.join(format!("run{index}.bam"));
    let mut writer = RecordWriter::create_with_workers(&path, Format::Bam, header, workers)?;
    for record in buffer.iter() {
        writer.write(header, record)?;
    }
    writer.finish()?;
    buffer.clear();
    Ok(path)
}

/// One entry per run in the merge heap. `Ord` is reversed so `BinaryHeap`, which is a
/// max-heap, yields the smallest name first.
struct Entry<T> {
    record: T,
    run: usize,
}

impl<T: Sortable> PartialEq for Entry<T> {
    fn eq(&self, other: &Self) -> bool {
        compare(&self.record, &other.record) == Ordering::Equal
    }
}
impl<T: Sortable> Eq for Entry<T> {}
impl<T: Sortable> PartialOrd for Entry<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl<T: Sortable> Ord for Entry<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed for the max-heap, then by run index so a tie resolves stably and the
        // output does not depend on how the input happened to be split into runs.
        compare(&other.record, &self.record).then_with(|| other.run.cmp(&self.run))
    }
}

fn merge(runs: &[PathBuf], dst: &Path, header: &Header, workers: NonZero<usize>) -> Result<()> {
    let mut readers: Vec<RecordReader> = runs
        .iter()
        .map(|p| RecordReader::open_with_workers(p, workers))
        .collect::<Result<_>>()?;

    let mut heap = BinaryHeap::new();
    for (run, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = reader.next().transpose()? {
            heap.push(Entry { record, run });
        }
    }

    let mut writer = RecordWriter::create_with_workers(dst, Format::Bam, header, workers)?;
    while let Some(Entry { record, run }) = heap.pop() {
        writer.write(header, &record)?;
        if let Some(next) = readers[run].next().transpose()? {
            heap.push(Entry { record: next, run });
        }
    }
    writer.finish()
}

/// What ordering needs from a record, and no more.
///
/// Implemented for both the decoded `RecordBuf` and the raw `bam::Record`, so the
/// comparator is written once and the fast path cannot drift from the slow one.
trait Sortable {
    fn sort_name(&self) -> &[u8];
    fn sort_flags(&self) -> noodles_sam::alignment::record::Flags;
}

impl Sortable for RecordBuf {
    fn sort_name(&self) -> &[u8] {
        self.name().map(|n| n.as_ref()).unwrap_or(&b""[..])
    }
    fn sort_flags(&self) -> noodles_sam::alignment::record::Flags {
        self.flags()
    }
}

impl Sortable for noodles_bam::Record {
    fn sort_name(&self) -> &[u8] {
        self.name().map(|n| n.as_ref()).unwrap_or(&b""[..])
    }
    fn sort_flags(&self) -> noodles_sam::alignment::record::Flags {
        self.flags()
    }
}

/// Order two records the way `samtools sort -n` would: by name, then by the READ1/READ2
/// flags so mates keep a stable relative order.
fn compare<T: Sortable>(a: &T, b: &T) -> Ordering {
    strnum_cmp(a.sort_name(), b.sort_name()).then_with(|| mate_rank(a).cmp(&mate_rank(b)))
}

fn mate_rank<T: Sortable>(record: &T) -> u8 {
    let flags = record.sort_flags();
    if flags.is_first_segment() {
        0
    } else if flags.is_last_segment() {
        1
    } else {
        2
    }
}

/// Sort a BAM by name without decoding its records.
///
/// The decoded path exists because the tagger needs decoded records. Sorting does not: it
/// needs a name and two flag bits, both of which `bam::Record` reads straight out of the
/// raw buffer it already holds. Skipping the decode-and-re-encode round trip is worth
/// roughly 3x, which is the difference between "slower than samtools" and "comparable to
/// it", so the BAM-to-BAM case gets its own path.
///
/// The comparator, the tie-breaking and the run/merge structure are shared with the decoded
/// path, so the two cannot order records differently.
fn sort_bam_raw(
    pool: &rayon::ThreadPool,
    src: &Path,
    dst: &Path,
    header: &Header,
    max_records_in_memory: usize,
    workers: NonZero<usize>,
) -> Result<()> {
    let spill_dir = tempfile::tempdir()?;
    let mut runs: Vec<PathBuf> = Vec::new();
    let mut buffer: Vec<noodles_bam::Record> = Vec::with_capacity(max_records_in_memory);

    let mut reader = open_raw(src, workers)?;
    let mut record = noodles_bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        // Move rather than clone: the reader refills a fresh buffer either way, and cloning
        // 2 million records is measurable.
        buffer.push(std::mem::take(&mut record));
        if buffer.len() >= max_records_in_memory {
            runs.push(spill_raw(
                pool,
                &mut buffer,
                header,
                spill_dir.path(),
                runs.len(),
                workers,
            )?);
        }
    }

    if runs.is_empty() {
        pool.install(|| buffer.par_sort_by(compare));
        let mut writer = create_raw(dst, header, workers)?;
        for record in &buffer {
            writer.write_record(header, record)?;
        }
        return finish_raw(writer);
    }

    if !buffer.is_empty() {
        runs.push(spill_raw(
            pool,
            &mut buffer,
            header,
            spill_dir.path(),
            runs.len(),
            workers,
        )?);
    }

    merge_raw(&runs, dst, header, workers)
}

type RawReader = noodles_bam::io::Reader<noodles_bgzf::io::MultithreadedReader<std::fs::File>>;
type RawWriter = noodles_bam::io::Writer<noodles_bgzf::io::MultithreadedWriter<std::fs::File>>;

fn open_raw(path: &Path, workers: NonZero<usize>) -> Result<RawReader> {
    use anyhow::Context;
    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open file '{}'", path.display()))?;
    let decoder = noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, file);
    let mut reader = noodles_bam::io::Reader::from(decoder);
    // The header has to be consumed before records, even when its contents are not wanted:
    // it is part of the stream, not a side channel.
    reader.read_header()?;
    Ok(reader)
}

fn create_raw(path: &Path, header: &Header, workers: NonZero<usize>) -> Result<RawWriter> {
    use anyhow::Context;
    let file = std::fs::File::create(path)
        .with_context(|| format!("Failed to write to file '{}'", path.display()))?;
    let encoder = noodles_bgzf::io::MultithreadedWriter::with_worker_count(workers, file);
    let mut writer = noodles_bam::io::Writer::from(encoder);
    writer.write_header(header)?;
    Ok(writer)
}

fn finish_raw(writer: RawWriter) -> Result<()> {
    use anyhow::Context;
    let mut encoder = writer.into_inner();
    encoder.finish().context("Failed to finalise BAM output")?;
    Ok(())
}

fn spill_raw(
    pool: &rayon::ThreadPool,
    buffer: &mut Vec<noodles_bam::Record>,
    header: &Header,
    dir: &Path,
    index: usize,
    workers: NonZero<usize>,
) -> Result<PathBuf> {
    pool.install(|| buffer.par_sort_by(compare));
    let path = dir.join(format!("run{index}.bam"));
    let mut writer = create_raw(&path, header, workers)?;
    for record in buffer.iter() {
        writer.write_record(header, record)?;
    }
    finish_raw(writer)?;
    buffer.clear();
    Ok(path)
}

fn merge_raw(runs: &[PathBuf], dst: &Path, header: &Header, workers: NonZero<usize>) -> Result<()> {
    let mut readers: Vec<RawReader> = runs
        .iter()
        .map(|p| open_raw(p, workers))
        .collect::<Result<_>>()?;

    let mut heap: BinaryHeap<Entry<noodles_bam::Record>> = BinaryHeap::new();
    for (run, reader) in readers.iter_mut().enumerate() {
        let mut record = noodles_bam::Record::default();
        if reader.read_record(&mut record)? != 0 {
            heap.push(Entry { record, run });
        }
    }

    let mut writer = create_raw(dst, header, workers)?;
    while let Some(Entry { record, run }) = heap.pop() {
        writer.write_record(header, &record)?;
        let mut next = noodles_bam::Record::default();
        if readers[run].read_record(&mut next)? != 0 {
            heap.push(Entry { record: next, run });
        }
    }
    finish_raw(writer)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The expected order is `samtools sort -n` (1.24) output for exactly these names, not a
    /// guess. If you change the comparator, re-derive this from samtools rather than
    /// adjusting it to match the code.
    const SAMTOOLS_ORDER: [&str; 8] = [
        "a1b2", "a1b10", "read", "read2", "read007", "read7", "read10", "readA",
    ];

    #[test]
    fn sorting_reproduces_the_samtools_order() {
        let mut names = SAMTOOLS_ORDER;
        names.reverse();
        names.sort_by(|x, y| strnum_cmp(x.as_bytes(), y.as_bytes()));
        assert_eq!(names, SAMTOOLS_ORDER);
    }

    #[test]
    fn digit_runs_compare_numerically_not_lexically() {
        assert_eq!(strnum_cmp(b"read2", b"read10"), Ordering::Less);
        assert_eq!(strnum_cmp(b"read10", b"read2"), Ordering::Greater);
        assert_eq!(strnum_cmp(b"read2", b"read2"), Ordering::Equal);
    }

    /// Equal value, different padding, is a tie broken by the padding rather than a tie.
    /// A comparator that returned Equal here would reorder mates non-deterministically.
    #[test]
    fn leading_zeroes_break_the_tie_rather_than_vanishing() {
        assert_eq!(strnum_cmp(b"read007", b"read7"), Ordering::Less);
        assert_eq!(strnum_cmp(b"read7", b"read007"), Ordering::Greater);
        assert_eq!(strnum_cmp(b"read007", b"read8"), Ordering::Less);
    }

    #[test]
    fn text_segments_compare_bytewise() {
        assert_eq!(strnum_cmp(b"readA", b"readB"), Ordering::Less);
        assert_eq!(strnum_cmp(b"HWI:1:2", b"HWI:1:10"), Ordering::Less);
    }

    #[test]
    fn a_prefix_sorts_before_what_extends_it() {
        assert_eq!(strnum_cmp(b"read", b"read1"), Ordering::Less);
        assert_eq!(strnum_cmp(b"", b"a"), Ordering::Less);
        assert_eq!(strnum_cmp(b"", b""), Ordering::Equal);
    }

    #[test]
    fn mixed_runs_alternate_correctly() {
        assert_eq!(strnum_cmp(b"a1b2", b"a1b10"), Ordering::Less);
        assert_eq!(strnum_cmp(b"a2b1", b"a10b1"), Ordering::Less);
    }
}
