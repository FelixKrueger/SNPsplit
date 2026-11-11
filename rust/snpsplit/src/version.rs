//! The three `--version` banners: the text Perl SNPsplit v0.9.2 prints, with `@VERSION@` in
//! place of its version number.
//!
//! They are `include_str!`d from `banners/` rather than written out here, so they cannot
//! drift from what the Perl actually prints.

use crate::Tool;

/// The suite version.
///
/// Taken from the package manifest rather than from a `VERSION` file beside it. The file was
/// the obvious "single source of truth" and turned out to be unpublishable: `include_str!`
/// cannot reach outside the package directory, so `cargo package` refused the crate. The
/// manifest is a source of truth that ships with the code.
pub const SUITE_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Stands for the version in each banner template; it occurs exactly once per template.
const PLACEHOLDER: &str = "@VERSION@";

const BANNER_SNPSPLIT: &str = include_str!("../banners/snpsplit.txt");
const BANNER_TAG2SORT: &str = include_str!("../banners/tag2sort.txt");
const BANNER_GENOME_PREP: &str = include_str!("../banners/genome_prep.txt");

fn template(tool: Tool) -> &'static str {
    match tool {
        Tool::Tag => BANNER_SNPSPLIT,
        Tool::Sort => BANNER_TAG2SORT,
        Tool::Prepare => BANNER_GENOME_PREP,
    }
}

/// The exact text the version flag writes to stdout for `tool`.
pub fn banner(tool: Tool) -> String {
    template(tool).replacen(PLACEHOLDER, SUITE_VERSION, 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOOLS: [Tool; 3] = [Tool::Tag, Tool::Sort, Tool::Prepare];

    #[test]
    fn every_template_has_the_placeholder_exactly_once() {
        for tool in TOOLS {
            assert_eq!(
                template(tool).matches(PLACEHOLDER).count(),
                1,
                "{tool:?} banner template"
            );
        }
    }

    /// The whole version token, not a substring.
    #[test]
    fn every_banner_states_exactly_the_suite_version() {
        for tool in TOOLS {
            let text = banner(tool);
            let line = text
                .lines()
                .find(|l| l.contains("ersion: "))
                .unwrap_or_else(|| panic!("{tool:?} banner has no version line"));
            let token = line.split("ersion: ").nth(1).unwrap().trim();
            assert_eq!(token, SUITE_VERSION, "{tool:?} banner");
        }
    }

    /// Perl prints these with leading and trailing blank lines that the fixtures diff.
    #[test]
    fn banners_keep_their_surrounding_whitespace() {
        assert!(banner(Tool::Tag).starts_with('\n'));
        assert!(banner(Tool::Sort).starts_with("\n\n"));
        assert!(banner(Tool::Prepare).starts_with('\t'));
        for tool in TOOLS {
            assert!(
                banner(tool).ends_with("\n\n"),
                "{tool:?} lost its trailing blank line"
            );
        }
    }
}
