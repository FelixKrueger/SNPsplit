//! Command line handling for `SNPsplit_genome_preparation`.
//!
//! The Perl uses Getopt::Long, whose diagnostics are part of the observable behaviour: a bad
//! option prints `Unknown option: name` and then `Please respecify command line options`.
//! That shape is reproduced rather than replaced with a nicer parser's output.

use std::path::{Path, PathBuf};

/// A run that ended before doing any work, and how it ended.
pub enum Stop {
    /// Printed something and exited successfully, like `--list_strains`.
    Ok,
    /// Died with this message, already terminated the way Perl's `die` would print it.
    Died(String),
}

/// Everything `process_commandline` resolves before the run proper begins.
#[derive(Debug, Default)]
pub struct Options {
    pub vcf_file: Option<PathBuf>,
    pub strain: Option<String>,
    pub strain2: Option<String>,
    pub list_strains: bool,
    pub skip_filtering: bool,
    pub full_sequence: bool,
    pub nmasking: bool,
    pub no_nmasking: bool,
    pub dual_hybrid: bool,
    pub genome_folder: Option<String>,
    pub genome_build: Option<String>,
    pub v7: bool,
    pub help: bool,
    pub version: bool,
    pub parallel: usize,
    pub download: bool,
    pub download_dir: Option<String>,
    pub ensembl_release: Option<u32>,
}

/// Parse the command line the way Getopt::Long does for this option set.
///
/// Returns `Err` with the Getopt::Long message on an unknown option or a missing value, which
/// the caller turns into the `Please respecify command line options` die.
pub fn parse(args: &[String]) -> Result<Options, String> {
    let mut opts = Options {
        parallel: crate::io::default_parallel(),
        ..Options::default()
    };
    let mut i = 0;

    while i < args.len() {
        let arg = &args[i];
        let Some(body) = arg.strip_prefix("--").or_else(|| arg.strip_prefix("-")) else {
            // Getopt::Long leaves non-options in @ARGV; this tool ignores them.
            i += 1;
            continue;
        };

        let (name, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v.to_string())),
            None => (body, None),
        };

        // A value-taking option consumes the next argument when there is no `=value`.
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
            "versions" => opts.version = true,
            "strain" => opts.strain = Some(take_value(&mut i)?),
            "strain2" => opts.strain2 = Some(take_value(&mut i)?),
            "list_strains" => opts.list_strains = true,
            "skip_filtering" => opts.skip_filtering = true,
            "vcf_file" => opts.vcf_file = Some(PathBuf::from(take_value(&mut i)?)),
            "full_sequence" => opts.full_sequence = true,
            "nmasking" => opts.nmasking = true,
            "no_nmasking" => opts.no_nmasking = true,
            "dual_hybrid" => opts.dual_hybrid = true,
            "reference_genome" => opts.genome_folder = Some(take_value(&mut i)?),
            "genome_build" => opts.genome_build = Some(take_value(&mut i)?),
            "v7_VCF" => opts.v7 = true,
            "parallel" => opts.parallel = take_value(&mut i)?.parse().unwrap_or(1),
            // The port's one addition, documented in rust/DESIGN.md. Nothing here contacts
            // the network unless --download is given.
            "download" => opts.download = true,
            "download_dir" => opts.download_dir = Some(take_value(&mut i)?),
            "ensembl_release" => {
                opts.ensembl_release = Some(
                    take_value(&mut i)?
                        .parse()
                        .map_err(|_| "Option ensembl_release requires a number".to_string())?,
                );
            }
            other => return Err(format!("Unknown option: {other}")),
        }

        i += 1;
    }

    Ok(opts)
}

/// The resolved configuration a run works from.
pub struct Config {
    pub vcf_file: Option<PathBuf>,
    pub strain: String,
    pub strain2: Option<String>,
    pub strain_index: Option<usize>,
    pub strain2_index: Option<usize>,
    pub genome_folder: String,
    pub skip_filtering: bool,
    pub nmasking: bool,
    pub full_sequence: bool,
    pub dual_hybrid: bool,
    pub genome_build: String,
    pub v7: bool,
    pub parallel: usize,
}

/// Whether a VCF path names the v7 combined SNP and INDEL release.
pub fn is_v7_release(path: &Path) -> bool {
    path.to_string_lossy()
        .contains("mgp_REL2005_snps_indels.vcf.gz")
}

/// The basename of a path, the way the Perl strips path information with `s/.*\///`.
pub fn basename(path: &Path) -> String {
    let s = path.to_string_lossy();
    match s.rfind('/') {
        Some(i) => s[i + 1..].to_string(),
        None => s.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_the_shape_the_fixtures_use() {
        let o = parse(&args(&[
            "--vcf_file",
            "snps.vcf",
            "--reference_genome",
            "genome",
            "--strain",
            "STRAIN_A",
        ]))
        .unwrap();
        assert_eq!(o.vcf_file.unwrap().to_str().unwrap(), "snps.vcf");
        assert_eq!(o.genome_folder.unwrap(), "genome");
        assert_eq!(o.strain.unwrap(), "STRAIN_A");
    }

    #[test]
    fn accepts_the_equals_form_as_getopt_long_does() {
        let o = parse(&args(&["--strain=STRAIN_A", "--genome_build=TESTBUILD"])).unwrap();
        assert_eq!(o.strain.unwrap(), "STRAIN_A");
        assert_eq!(o.genome_build.unwrap(), "TESTBUILD");
    }

    /// A value containing a space arrives as one argv entry, so nothing has to quote it. The
    /// Perl reaches gunzip through a shell-free open for the same reason.
    #[test]
    fn a_value_may_contain_spaces() {
        let o = parse(&args(&["--genome_build", "test build"])).unwrap();
        assert_eq!(o.genome_build.unwrap(), "test build");
    }

    #[test]
    fn the_download_options_parse() {
        let o = parse(&args(&[
            "--download",
            "--download_dir",
            "refs",
            "--ensembl_release",
            "115",
        ]))
        .unwrap();
        assert!(o.download);
        assert_eq!(o.download_dir.unwrap(), "refs");
        assert_eq!(o.ensembl_release.unwrap(), 115);
    }

    /// Without the flag the tool has no reason to touch the network, and the other two
    /// options are inert.
    #[test]
    fn downloading_is_off_unless_asked_for() {
        let o = parse(&args(&["--download_dir", "refs"])).unwrap();
        assert!(!o.download);
    }

    #[test]
    fn an_unknown_option_reports_it_the_getopt_long_way() {
        let err = parse(&args(&["--frobnicate"])).unwrap_err();
        assert_eq!(err, "Unknown option: frobnicate");
    }

    #[test]
    fn the_v7_release_is_recognised_by_name_anywhere_in_the_path() {
        assert!(is_v7_release(Path::new("mgp_REL2005_snps_indels.vcf.gz")));
        assert!(is_v7_release(Path::new(
            "/data/vcf dir/mgp_REL2005_snps_indels.vcf.gz"
        )));
        assert!(!is_v7_release(Path::new("mgp_REL2021_snps.vcf.gz")));
    }

    #[test]
    fn basename_strips_directories_the_way_the_perl_regex_does() {
        assert_eq!(basename(Path::new("/a/b/c.vcf")), "c.vcf");
        assert_eq!(basename(Path::new("c.vcf")), "c.vcf");
        assert_eq!(basename(Path::new("/a/vcf dir/c.vcf.gz")), "c.vcf.gz");
    }
}
