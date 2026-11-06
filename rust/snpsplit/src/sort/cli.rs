//! Command line handling for `tag2sort`.

use std::path::PathBuf;

/// Everything `process_commandline` resolves.
#[derive(Debug, Default)]
pub struct Options {
    pub paired: bool,
    pub hic: bool,
    pub verbose: bool,
    pub samtools_path: Option<String>,
    pub sam: bool,
    pub conflict: bool,
    pub singletons: bool,
    pub output_dir: Option<String>,
    pub parent_dir: Option<String>,
    pub help: bool,
    pub version: bool,
    pub parallel: usize,
    pub files: Vec<PathBuf>,
}

/// Parse the command line the way Getopt::Long does for this option set.
pub fn parse(args: &[String]) -> Result<Options, String> {
    let mut opts = Options::default();
    let mut i = 0;

    while i < args.len() {
        let arg = &args[i];
        let Some(body) = arg.strip_prefix("--").or_else(|| arg.strip_prefix("-")) else {
            opts.files.push(PathBuf::from(arg));
            i += 1;
            continue;
        };

        let (name, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v.to_string())),
            None => (body, None),
        };

        let take_value = |i: &mut usize| -> Result<String, String> {
            if let Some(v) = inline.clone() {
                return Ok(v);
            }
            *i += 1;
            args.get(*i)
                .cloned()
                .ok_or_else(|| format!("Option {name} requires an argument"))
        };

        match name {
            "help" | "man" => opts.help = true,
            // Singular, and the plural is rejected. See src/tool.rs.
            "version" => opts.version = true,
            "output_dir" => opts.output_dir = Some(take_value(&mut i)?),
            "paired" => opts.paired = true,
            "verbose" => opts.verbose = true,
            "hic" => opts.hic = true,
            "samtools_path" => opts.samtools_path = Some(take_value(&mut i)?),
            "sam" => opts.sam = true,
            "conflicting" | "weird" => opts.conflict = true,
            "singletons" => opts.singletons = true,
            "parent_dir" => opts.parent_dir = Some(take_value(&mut i)?),
            "parallel" => {
                opts.parallel = take_value(&mut i)?.parse().unwrap_or(1);
            }
            other => return Err(format!("Unknown option: {other}")),
        }

        i += 1;
    }

    Ok(opts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn collects_input_files_and_flags() {
        let o = parse(&args(&[
            "--paired",
            "--conflicting",
            "sample.allele_flagged.bam",
        ]))
        .unwrap();
        assert!(o.paired);
        assert!(o.conflict);
        assert_eq!(o.files.len(), 1);
    }

    /// The Perl spells this option two ways and both are in the wild.
    #[test]
    fn conflicting_and_weird_are_the_same_option() {
        assert!(parse(&args(&["--weird"])).unwrap().conflict);
        assert!(parse(&args(&["--conflicting"])).unwrap().conflict);
    }

    #[test]
    fn an_empty_output_dir_is_kept_as_empty_rather_than_dropped() {
        // SNPsplit passes --output_dir '' when it has no directory to add, and the Perl
        // deliberately leaves that alone rather than turning it into '/'.
        let o = parse(&args(&["--output_dir", ""])).unwrap();
        assert_eq!(o.output_dir.as_deref(), Some(""));
    }
}
