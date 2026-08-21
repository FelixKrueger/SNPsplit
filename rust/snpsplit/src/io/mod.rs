//! Alignment I/O. Pure Rust via noodles: no samtools subprocess, no htslib.

mod namesort;
mod reader;
mod writer;

use std::num::NonZero;

pub use namesort::{sort_by_name, sort_by_name_with_workers, strnum_cmp};
pub use reader::{Format, RecordReader, sniff};
pub use writer::{RecordWriter, add_pg_line, command_line, render_sam_line};

/// Resolve a `--parallel N` value into a worker count.
///
/// `0` means every available core, which is the only value whose meaning is not simply the
/// number itself. Anything else is taken literally, including counts above the core count:
/// oversubscription is the caller's business, and a scheduler that hands out four cores does
/// not want us second-guessing it.
pub fn worker_count(requested: usize) -> NonZero<usize> {
    const ONE: NonZero<usize> = NonZero::new(1).unwrap();

    match NonZero::new(requested) {
        Some(n) => n,
        None => std::thread::available_parallelism().unwrap_or(ONE),
    }
}

/// The default worker count when `--parallel` was not given.
///
/// `SNPSPLIT_PARALLEL` sets it, which is how a container or a scheduler can hand every tool
/// in a pipeline the same core budget without rewriting the command lines. The flag always
/// wins when it is present, and an unparseable value is ignored rather than fatal: a stray
/// environment variable should not stop a run.
///
/// It is also what lets the fixture suites run at more than one worker count without the
/// runners growing a way to pass extra arguments.
pub fn default_parallel() -> usize {
    std::env::var("SNPSPLIT_PARALLEL")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_literal_count_is_taken_literally() {
        assert_eq!(worker_count(1).get(), 1);
        assert_eq!(worker_count(4).get(), 4);
        // Oversubscription is allowed: the caller may know something we do not.
        assert_eq!(worker_count(1024).get(), 1024);
    }

    #[test]
    fn zero_means_every_core() {
        let all = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        assert_eq!(worker_count(0).get(), all);
    }
}
