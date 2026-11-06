//! Alignment I/O. Pure Rust via noodles: no samtools subprocess, no htslib.

mod namesort;
mod reader;
mod writer;

use std::num::NonZero;

pub use namesort::{sort_by_name, sort_by_name_with_workers, strnum_cmp};
pub use reader::{Format, RecordReader, sniff};
pub use writer::{RecordWriter, add_pg_line, command_line};

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
