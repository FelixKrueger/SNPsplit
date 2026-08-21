//! Resolving long option names the way Getopt::Long does by default.
//!
//! Two of its defaults are load-bearing and neither is obvious: `ignore_case` makes
//! `--snp_file` and `--SNP_file` the same option, and `auto_abbrev` accepts any unambiguous
//! prefix, so `--snp` works too. Both are reachable from any existing script, and the first
//! one matters more than it sounds: SNPsplit's own help documents `--snp_file` while its
//! code registers `SNP_file`, so a parser that only accepted the registered spelling would
//! reject the spelling the tool tells people to use.

/// Resolve `input` against the registered option names.
///
/// An exact match wins over a prefix, which is why `--strain` selects `strain` rather than
/// being ambiguous with `strain2`.
pub fn resolve<'a>(input: &str, names: &[&'a str]) -> Result<&'a str, String> {
    resolve_with(input, names, &[])
}

/// Resolve `input`, preferring an option the Perl already had over one this port added.
///
/// Adding an option can make an abbreviation ambiguous that was not: `--p` meant `paired`
/// until `--parallel` existed. Reporting that as ambiguous would break a working script for
/// no reason the user can see, and silently picking the new option would be worse, so a tie
/// between an original option and an added one goes to the original. `--par` still reaches
/// `parallel`, because that is not a tie.
///
/// A tie between two original options is ambiguous exactly as it always was.
pub fn resolve_with<'a>(input: &str, names: &[&'a str], added: &[&str]) -> Result<&'a str, String> {
    let lowered = input.to_ascii_lowercase();

    if let Some(name) = names
        .iter()
        .find(|name| name.to_ascii_lowercase() == lowered)
    {
        return Ok(name);
    }

    let matches: Vec<&&str> = names
        .iter()
        .filter(|name| name.to_ascii_lowercase().starts_with(&lowered))
        .collect();

    if matches.len() > 1 {
        let original: Vec<&&&str> = matches
            .iter()
            .filter(|name| !added.contains(**name))
            .collect();
        if original.len() == 1 {
            return Ok(original[0]);
        }
    }

    match matches.len() {
        1 => Ok(matches[0]),
        0 => Err(format!("Unknown option: {input}")),
        _ => {
            // Getopt::Long lists the candidates lowercased and sorted.
            let mut candidates: Vec<String> = matches
                .iter()
                .map(|name| name.to_ascii_lowercase())
                .collect();
            candidates.sort();
            Err(format!(
                "Option {input} is ambiguous ({})",
                candidates.join(", ")
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAGGER: &[&str] = &[
        "help",
        "man",
        "versions",
        "output_dir",
        "SNP_file",
        "no_sorting",
        "verbose",
        "samtools_path",
        "sam",
        "paired",
        "single_end",
        "hic",
        "conflicting",
        "weird",
        "singletons",
        "bisulfite",
        "skip_tag2sort",
    ];

    #[test]
    fn an_exact_name_resolves_to_itself() {
        assert_eq!(resolve("SNP_file", TAGGER).unwrap(), "SNP_file");
        assert_eq!(resolve("paired", TAGGER).unwrap(), "paired");
    }

    /// The spelling SNPsplit's own help documents. A case-sensitive parser rejects it.
    #[test]
    fn the_documented_lowercase_spelling_resolves() {
        assert_eq!(resolve("snp_file", TAGGER).unwrap(), "SNP_file");
        assert_eq!(resolve("SNP_FILE", TAGGER).unwrap(), "SNP_file");
    }

    #[test]
    fn an_unambiguous_prefix_resolves() {
        assert_eq!(resolve("snp", TAGGER).unwrap(), "SNP_file");
        assert_eq!(resolve("bisul", TAGGER).unwrap(), "bisulfite");
    }

    /// An exact name that is also a prefix of another is not ambiguous.
    #[test]
    fn an_exact_name_wins_over_being_a_prefix() {
        let names = &["strain", "strain2"];
        assert_eq!(resolve("strain", names).unwrap(), "strain");
        assert_eq!(resolve("strain2", names).unwrap(), "strain2");
    }

    /// The message and the ordering are Getopt::Long's, verified against it.
    #[test]
    fn an_ambiguous_prefix_lists_the_candidates_lowercased_and_sorted() {
        let err = resolve("s", TAGGER).unwrap_err();
        assert_eq!(
            err,
            "Option s is ambiguous (sam, samtools_path, single_end, singletons, skip_tag2sort, snp_file)"
        );

        let err = resolve("str", &["strain", "strain2"]).unwrap_err();
        assert_eq!(err, "Option str is ambiguous (strain, strain2)");
    }

    /// Adding an option must not change what an abbreviation meant before it existed.
    #[test]
    fn an_added_option_does_not_steal_an_existing_abbreviation() {
        let names = &["paired", "parallel"];
        let added = &["parallel"];
        assert_eq!(resolve_with("p", names, added).unwrap(), "paired");
        // Not a tie, so the added option is still reachable by prefix.
        assert_eq!(resolve_with("par", names, added).unwrap(), "parallel");
        assert_eq!(resolve_with("parallel", names, added).unwrap(), "parallel");
    }

    /// A tie between two options the Perl already had is ambiguous exactly as it was.
    #[test]
    fn a_tie_between_original_options_is_still_ambiguous() {
        let err = resolve_with("strain", &["strainA", "strainB"], &["parallel"]).unwrap_err();
        assert_eq!(err, "Option strain is ambiguous (straina, strainb)");
    }

    #[test]
    fn an_unknown_option_is_reported_as_unknown() {
        assert_eq!(
            resolve("frobnicate", TAGGER).unwrap_err(),
            "Unknown option: frobnicate"
        );
    }
}
