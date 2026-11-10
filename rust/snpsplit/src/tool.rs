//! Resolving which of the three tools a given invocation means.

use crate::Tool;

/// The tool a classic-name invocation selects, if any.
///
/// The path is stripped first: pipelines invoke these by absolute path, and the fixture
/// runner stages copies into a scratch directory before running them.
pub fn from_argv0(argv0: &str) -> Option<Tool> {
    let name = argv0.rsplit('/').next().unwrap_or(argv0);
    match name {
        "SNPsplit" => Some(Tool::Tag),
        "tag2sort" => Some(Tool::Sort),
        "SNPsplit_genome_preparation" => Some(Tool::Prepare),
        _ => None,
    }
}

/// The tool a `snpsplit <subcommand>` invocation selects, if any.
pub fn from_subcommand(name: &str) -> Option<Tool> {
    match name {
        "tag" => Some(Tool::Tag),
        "sort" => Some(Tool::Sort),
        "prepare" => Some(Tool::Prepare),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classic_names_resolve() {
        assert_eq!(from_argv0("SNPsplit"), Some(Tool::Tag));
        assert_eq!(from_argv0("tag2sort"), Some(Tool::Sort));
        assert_eq!(
            from_argv0("SNPsplit_genome_preparation"),
            Some(Tool::Prepare)
        );
    }

    #[test]
    fn a_classic_name_resolves_through_a_path() {
        assert_eq!(from_argv0("/tmp/scratch/bin/SNPsplit"), Some(Tool::Tag));
    }

    #[test]
    fn the_multicall_name_itself_is_not_a_tool() {
        assert_eq!(from_argv0("snpsplit"), None);
        assert_eq!(from_argv0("/usr/local/bin/snpsplit"), None);
    }

    #[test]
    fn subcommands_resolve() {
        assert_eq!(from_subcommand("tag"), Some(Tool::Tag));
        assert_eq!(from_subcommand("sort"), Some(Tool::Sort));
        assert_eq!(from_subcommand("prepare"), Some(Tool::Prepare));
        assert_eq!(from_subcommand("split"), None);
    }
}
