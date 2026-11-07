//! The three `--version` banners, captured byte for byte from Perl SNPsplit v0.9.0.
//!
//! They are `include_str!`d from `banners/` rather than written out here, so they cannot
//! drift from what the Perl actually prints. Regenerate them, do not edit them.

use crate::Tool;

/// The suite version, from `rust/VERSION`.
pub const SUITE_VERSION: &str = include_str!("../../VERSION").trim_ascii_end();

const BANNER_SNPSPLIT: &str = include_str!("../banners/snpsplit.txt");
const BANNER_TAG2SORT: &str = include_str!("../banners/tag2sort.txt");
const BANNER_GENOME_PREP: &str = include_str!("../banners/genome_prep.txt");

/// The exact text the version flag writes to stdout for `tool`.
pub fn banner(tool: Tool) -> &'static str {
    match tool {
        Tool::Tag => BANNER_SNPSPLIT,
        Tool::Sort => BANNER_TAG2SORT,
        Tool::Prepare => BANNER_GENOME_PREP,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The banners carry the version as literal text. If VERSION moves and the banners are
    /// not regenerated, the tool would report two different versions depending on the flag.
    #[test]
    fn every_banner_states_the_suite_version() {
        for tool in [Tool::Tag, Tool::Sort, Tool::Prepare] {
            assert!(
                banner(tool).contains(SUITE_VERSION),
                "{tool:?} banner does not mention version {SUITE_VERSION}; regenerate banners/",
            );
        }
    }

    /// Perl prints these with leading and trailing blank lines that the fixtures diff.
    #[test]
    fn banners_keep_their_surrounding_whitespace() {
        assert!(banner(Tool::Tag).starts_with('\n'));
        assert!(banner(Tool::Sort).starts_with("\n\n"));
        assert!(banner(Tool::Prepare).starts_with('\t'));
        for tool in [Tool::Tag, Tool::Sort, Tool::Prepare] {
            assert!(
                banner(tool).ends_with("\n\n"),
                "{tool:?} lost its trailing blank line"
            );
        }
    }
}
