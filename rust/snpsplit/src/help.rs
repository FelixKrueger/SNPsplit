//! The help pages: Perl SNPsplit v0.9.2's, plus each `ADDED` option and Perl options it left out.
//! A wrapped line never starts with `-`, since that is how a listed option is recognised.

use crate::Tool;

const HELP_SNPSPLIT: &str = include_str!("../help/snpsplit.txt");
const HELP_TAG2SORT: &str = include_str!("../help/tag2sort.txt");
const HELP_GENOME_PREP: &str = include_str!("../help/genome_prep.txt");

/// The help page for `tool`.
pub fn page(tool: Tool) -> &'static str {
    match tool {
        Tool::Tag => HELP_SNPSPLIT,
        Tool::Sort => HELP_TAG2SORT,
        Tool::Prepare => HELP_GENOME_PREP,
    }
}

/// The exit status `--help` produces.
///
/// Not the same for all three: the tagger and the sorter exit 1, the genome preparation exits
/// 0. That is an asymmetry in the Perl rather than a decision, and it is reproduced because a
/// script that tests `--help` for success would notice.
pub fn exit_status(tool: Tool) -> u8 {
    match tool {
        Tool::Tag | Tool::Sort => 1,
        Tool::Prepare => 0,
    }
}

/// Every `--option` a help page lists as an option of its own tool.
///
/// Only the option column counts: the first field of a line starting with `-`, split on `/`
/// (`-o/--output_dir`, `--conflicting/--weird`). Help text also mentions options in prose,
/// including options belonging to the other tools in the suite, and those are not promises about
/// what this tool accepts: tag2sort's page explains that `--hic` "sets the flags --paired and
/// --no_sort", neither of which tag2sort has.
pub fn documented_options(page: &str) -> Vec<String> {
    let mut found: Vec<String> = Vec::new();

    for line in page.lines() {
        let line = line.trim_start();
        if !line.starts_with('-') {
            continue;
        }
        let field = line.split_whitespace().next().unwrap_or("");
        for part in field.split('/') {
            let Some(rest) = part.strip_prefix("--") else {
                continue;
            };
            let name: String = rest
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                .collect();
            if !name.is_empty() && !found.contains(&name) {
                found.push(name);
            }
        }
    }

    found
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_page_is_present_and_names_its_tool() {
        assert!(page(Tool::Tag).contains("SNPsplit"));
        assert!(page(Tool::Sort).contains("tag2sort"));
        assert!(page(Tool::Prepare).contains("SNPsplit_genome_preparation"));
    }

    /// Every option a help page documents has to be one the tool accepts.
    ///
    /// This is the test that matters, and it started as a weaker one that asserted the help
    /// mentioned `--SNP_file`. It does not: it documents `--snp_file`, and the parser only
    /// took the registered spelling, so the tool rejected the spelling its own help told
    /// people to use. Checking the pages against the parsers is what catches that.
    #[test]
    fn every_documented_option_is_accepted_by_its_tool() {
        let cases: [(Tool, &[&str]); 3] = [
            (Tool::Tag, crate::tag::cli::NAMES),
            (Tool::Sort, crate::sort::cli::NAMES),
            (Tool::Prepare, crate::genome_prep::cli::NAMES),
        ];

        for (tool, names) in cases {
            for option in documented_options(page(tool)) {
                assert!(
                    crate::optmatch::resolve(&option, names).is_ok(),
                    "{tool:?} documents --{option} but does not accept it"
                );
            }
        }
    }

    /// Every option the port adds is listed on its tool's help page and accepted by its parser.
    #[test]
    fn every_added_option_is_listed_and_accepted() {
        let cases: [(Tool, &[&str], &[&str]); 3] = [
            (Tool::Tag, crate::tag::cli::ADDED, crate::tag::cli::NAMES),
            (Tool::Sort, crate::sort::cli::ADDED, crate::sort::cli::NAMES),
            (
                Tool::Prepare,
                crate::genome_prep::cli::ADDED,
                crate::genome_prep::cli::NAMES,
            ),
        ];

        for (tool, added, names) in cases {
            let listed = documented_options(page(tool));
            for name in added {
                assert!(
                    listed.iter().any(|o| o == name),
                    "{tool:?} adds --{name} but its help page does not list it"
                );
                assert!(
                    names.contains(name),
                    "{tool:?} adds --{name} but does not accept it"
                );
            }
        }
    }

    /// Every option a tool accepts is on its help page, or `UNLISTED` and not on the page.
    #[test]
    fn every_accepted_option_is_listed_or_unlisted() {
        // Tool, NAMES, ADDED, UNLISTED.
        type Case = (
            Tool,
            &'static [&'static str],
            &'static [&'static str],
            &'static [&'static str],
        );
        let cases: [Case; 3] = [
            (
                Tool::Tag,
                crate::tag::cli::NAMES,
                crate::tag::cli::ADDED,
                crate::tag::cli::UNLISTED,
            ),
            (
                Tool::Sort,
                crate::sort::cli::NAMES,
                crate::sort::cli::ADDED,
                crate::sort::cli::UNLISTED,
            ),
            (
                Tool::Prepare,
                crate::genome_prep::cli::NAMES,
                crate::genome_prep::cli::ADDED,
                crate::genome_prep::cli::UNLISTED,
            ),
        ];

        for (tool, names, added, unlisted) in cases {
            let listed: Vec<&str> = documented_options(page(tool))
                .iter()
                .map(|o| {
                    crate::optmatch::resolve_with(o, names, added)
                        .unwrap_or_else(|e| panic!("{tool:?} documents --{o}: {e}"))
                })
                .collect();
            for name in names {
                assert!(
                    listed.contains(name) || unlisted.contains(name),
                    "{tool:?} accepts --{name}, which is neither on its help page nor UNLISTED"
                );
            }
            for name in unlisted {
                assert!(
                    names.contains(name),
                    "{tool:?} lists --{name} as UNLISTED but does not accept it"
                );
                assert!(
                    !listed.contains(name),
                    "{tool:?} lists --{name} as UNLISTED but its page lists it"
                );
            }
        }
    }

    #[test]
    fn option_extraction_takes_the_list_and_not_the_prose() {
        let found =
            documented_options("    --snp_file <file>   text\n    --paired            more\n");
        assert_eq!(found, vec!["snp_file", "paired"]);

        let found = documented_options(
            "-o/--output_dir <dir>  text\n--conflicting/--weird  text\n\t--skip_tag2sort  text\n",
        );
        assert_eq!(
            found,
            vec!["output_dir", "conflicting", "weird", "skip_tag2sort"]
        );
        assert!(documented_options("- a bullet\n-- x\n--\n").is_empty());

        // Options named mid-sentence belong to whatever the sentence is about, which is not
        // necessarily this tool.
        assert!(documented_options("this also sets --paired and --no_sort").is_empty());
        assert!(documented_options("see -- below").is_empty());
    }

    #[test]
    fn the_exit_status_asymmetry_is_reproduced() {
        assert_eq!(exit_status(Tool::Tag), 1);
        assert_eq!(exit_status(Tool::Sort), 1);
        assert_eq!(exit_status(Tool::Prepare), 0);
    }
}
