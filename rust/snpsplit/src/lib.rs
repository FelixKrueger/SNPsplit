pub mod tool;
pub mod version;

/// The three tools the one binary provides.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    /// `SNPsplit`: tag alignments by allele.
    Tag,
    /// `tag2sort`: sort tagged alignments into per-genome files.
    Sort,
    /// `SNPsplit_genome_preparation`: build N-masked or full-sequence genomes.
    Prepare,
}
