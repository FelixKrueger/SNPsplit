use std::process::ExitCode;

use snpsplit::{Tool, genome_prep, sort, tag, tool, version};

fn main() -> ExitCode {
    let mut args: Vec<String> = std::env::args().collect();
    let argv0 = args.first().cloned().unwrap_or_default();

    // A classic name selects the tool and leaves the rest of the args alone. Otherwise the
    // first argument is a subcommand and is consumed.
    let selected = match tool::from_argv0(&argv0) {
        Some(t) => Some(t),
        None => match args.get(1).map(String::as_str) {
            Some("--version") => {
                println!("snpsplit {}", version::SUITE_VERSION);
                return ExitCode::SUCCESS;
            }
            Some("--help" | "-h") => {
                println!("{USAGE}");
                return ExitCode::SUCCESS;
            }
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

/// Run the selected tool.
fn run(tool: Tool, args: &[String]) -> ExitCode {
    match tool {
        Tool::Prepare => genome_prep::run(args),
        Tool::Sort => sort::run(args),
        Tool::Tag => tag::run(args),
    }
}

const USAGE: &str = "Usage: snpsplit <tag|sort|prepare> ...\n   \
                     or: invoke as SNPsplit, tag2sort, or SNPsplit_genome_preparation\n\
                     Each tool's own options: snpsplit <tag|sort|prepare> --help";

fn usage_error(argv0: &str, subcommand: Option<&str>) -> ExitCode {
    match subcommand {
        Some(sub) => eprintln!("snpsplit: unknown subcommand '{sub}'"),
        None => eprintln!("snpsplit: cannot tell which tool '{argv0}' means"),
    }
    eprintln!("{USAGE}");
    ExitCode::FAILURE
}
