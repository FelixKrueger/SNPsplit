//! Alignment I/O. Pure Rust via noodles: no samtools subprocess, no htslib.

mod reader;
mod writer;

pub use reader::{Format, RecordReader, sniff};
pub use writer::RecordWriter;
