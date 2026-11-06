use std::process::ExitCode;

use snpsplit::{Tool, tool, version};

fn main() -> ExitCode {
    let mut args: Vec<String> = std::env::args().collect();
    let argv0 = args.first().cloned().unwrap_or_default();

    // A classic name selects the tool and leaves the rest of the args alone. Otherwise the
    // first argument is a subcommand and is consumed.
    let selected = match tool::from_argv0(&argv0) {
        Some(t) => Some(t),
        None => match args.get(1) {
            Some(sub) => {
                let t = tool::from_subcommand(sub);
                if t.is_some() {
                    args.remove(1);
                }
                t
            }
            None => None,
        },
    };

    let Some(selected) = selected else {
        return usage_error(&argv0, args.get(1).map(String::as_str));
    };

    run(selected, &args[1..])
}

/// Phase 0 stops here: the tools themselves land in later PRs. Everything except the
/// version flag exits non-zero with a message that names the tool, so a fixture that
/// reaches an unported path fails loudly rather than producing an empty output tree.
fn run(tool: Tool, args: &[String]) -> ExitCode {
    let flag = format!("--{}", tool::version_flag(tool));

    if args.iter().any(|a| a == &flag) {
        print!("{}", version::banner(tool));
        return ExitCode::SUCCESS;
    }

    // The Perl Getopt::Long rejection, reproduced: an unknown long option is fatal, and the
    // near-miss spelling of the version flag is the one users actually hit.
    for arg in args {
        if let Some(name) = arg.strip_prefix("--")
            && (name == "version" || name == "versions")
        {
            eprintln!("Unknown option: {name}");
            eprintln!("Please respecify command line options");
            return ExitCode::from(255);
        }
    }

    eprintln!("{tool:?} is not implemented in this build yet");
    ExitCode::FAILURE
}

fn usage_error(argv0: &str, subcommand: Option<&str>) -> ExitCode {
    match subcommand {
        Some(sub) => eprintln!("snpsplit: unknown subcommand '{sub}'"),
        None => eprintln!("snpsplit: cannot tell which tool '{argv0}' means"),
    }
    eprintln!("Usage: snpsplit <tag|sort|prepare> ...");
    eprintln!("   or: invoke as SNPsplit, tag2sort, or SNPsplit_genome_preparation");
    ExitCode::FAILURE
}
