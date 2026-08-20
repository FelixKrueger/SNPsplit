//! Alignment I/O. Pure Rust via noodles: no samtools subprocess, no htslib.

mod namesort;
mod reader;
mod writer;

pub use namesort::{sort_by_name, strnum_cmp};
pub use reader::{Format, RecordReader, sniff};
pub use writer::RecordWriter;
