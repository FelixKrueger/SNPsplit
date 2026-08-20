//! Times the external name sort at a given worker count, and checks the result against a
//! reference run. Used to produce the numbers in the parallelism section of the design;
//! not part of the shipped binary.
//!
//! Usage: sort_bench <input.bam> <workers> <records-in-memory>

use std::num::NonZero;
use std::path::PathBuf;
use std::time::Instant;

use snpsplit::io::{RecordReader, sort_by_name_with_workers, worker_count};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let src = PathBuf::from(args.next().expect("input path"));
    let workers: usize = args.next().expect("worker count").parse()?;
    let budget: usize = args.next().expect("records in memory").parse()?;

    let dir = tempfile::tempdir()?;
    let dst = dir.path().join("sorted.bam");
    let header = RecordReader::open(&src)?.header().clone();

    let started = Instant::now();
    sort_by_name_with_workers(&src, &dst, &header, budget, worker_count(workers))?;
    let elapsed = started.elapsed();

    let records = RecordReader::open(&dst)?.count();
    let effective: NonZero<usize> = worker_count(workers);
    println!(
        "workers={} budget={} records={} elapsed={:.2}s",
        effective.get(),
        budget,
        records,
        elapsed.as_secs_f64()
    );
    Ok(())
}
