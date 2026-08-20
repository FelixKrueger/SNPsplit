//! Alignment I/O. Pure Rust via noodles: no samtools subprocess, no htslib.

mod reader;

pub use reader::{Format, RecordReader, sniff};
