//! Multicall behaviour, exercised through the real binary.
//!
//! `argv[0]` dispatch cannot be unit-tested: it depends on how the process was actually
//! spawned. These run the built binary under each of its names.

use std::os::unix::fs::symlink;
use std::path::{Path, PathBuf};
use std::process::Command;

use assert_cmd::cargo::cargo_bin;
use tempfile::TempDir;

/// Links the built binary under `name` in a fresh directory and returns the link path.
fn linked_as(dir: &Path, name: &str) -> PathBuf {
    let link = dir.join(name);
    symlink(cargo_bin("snpsplit"), &link).expect("symlink the binary under its classic name");
    link
}

#[test]
fn snpsplit_prints_its_banner_for_versions() {
    let dir = TempDir::new().unwrap();
    let out = Command::new(linked_as(dir.path(), "SNPsplit"))
        .arg("--versions")
        .output()
        .unwrap();

    assert!(
        out.status.success(),
        "exit status was {:?}",
        out.status.code()
    );
    assert_eq!(
        String::from_utf8_lossy(&out.stdout),
        include_str!("../banners/snpsplit.txt"),
    );
}

#[test]
fn tag2sort_prints_its_banner_for_the_singular_version() {
    let dir = TempDir::new().unwrap();
    let out = Command::new(linked_as(dir.path(), "tag2sort"))
        .arg("--version")
        .output()
        .unwrap();

    assert!(out.status.success());
    assert_eq!(
        String::from_utf8_lossy(&out.stdout),
        include_str!("../banners/tag2sort.txt"),
    );
}

/// Perl: `Unknown option: versions` then `Please respecify command line options`, status 255.
/// Accepting the plural here would be a kindness that silently diverges from the Perl.
#[test]
fn tag2sort_rejects_the_plural_version_flag() {
    let dir = TempDir::new().unwrap();
    let out = Command::new(linked_as(dir.path(), "tag2sort"))
        .arg("--versions")
        .output()
        .unwrap();

    assert_eq!(out.status.code(), Some(255));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        stderr,
        "Unknown option: versions\nPlease respecify command line options\n"
    );
}

#[test]
fn genome_preparation_prints_its_banner() {
    let dir = TempDir::new().unwrap();
    let out = Command::new(linked_as(dir.path(), "SNPsplit_genome_preparation"))
        .arg("--versions")
        .output()
        .unwrap();

    assert!(out.status.success());
    assert_eq!(
        String::from_utf8_lossy(&out.stdout),
        include_str!("../banners/genome_prep.txt"),
    );
}

#[test]
fn the_subcommand_form_reaches_the_same_tool() {
    let out = Command::new(cargo_bin("snpsplit"))
        .args(["tag", "--versions"])
        .output()
        .unwrap();

    assert!(out.status.success());
    assert_eq!(
        String::from_utf8_lossy(&out.stdout),
        include_str!("../banners/snpsplit.txt"),
    );
}

#[test]
fn an_unknown_name_and_an_unknown_subcommand_both_fail_loudly() {
    let dir = TempDir::new().unwrap();
    let out = Command::new(linked_as(dir.path(), "SNPsplit_typo"))
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("SNPsplit_typo"));

    let out = Command::new(cargo_bin("snpsplit"))
        .arg("frobnicate")
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("frobnicate"));
}

/// The fixture runners copy the implementation into a scratch directory rather than
/// invoking it in place. Dispatch must therefore survive being copied under a classic name,
/// not merely being linked as one.
#[test]
fn dispatch_survives_being_copied_under_a_classic_name() {
    let dir = TempDir::new().unwrap();
    let copied = dir.path().join("SNPsplit");
    std::fs::copy(cargo_bin("snpsplit"), &copied).unwrap();

    let out = Command::new(&copied).arg("--versions").output().unwrap();

    assert!(out.status.success());
    assert_eq!(
        String::from_utf8_lossy(&out.stdout),
        include_str!("../banners/snpsplit.txt"),
    );
}

/// Every spelling of the version flag gives the Perl's exit status, and the Perl's banner when it succeeds.
#[test]
fn version_flags_behave_as_the_perl_does() {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let dir = TempDir::new().unwrap();
    for tool in ["SNPsplit", "tag2sort", "SNPsplit_genome_preparation"] {
        let rust = linked_as(dir.path(), tool);
        for flags in [
            &["--version"][..],
            &["--versions"],
            &["--vers"],
            &["--v"],
            &["--versions", "--bogus"],
            &["--version", "--bogus"],
        ] {
            let perl = Command::new("perl")
                .arg(repo.join(tool))
                .args(flags)
                .output()
                .unwrap();
            let ours = Command::new(&rust).args(flags).output().unwrap();
            assert_eq!(ours.status.code(), perl.status.code(), "{tool} {flags:?}");
            if perl.status.success() {
                assert_eq!(ours.stdout, perl.stdout, "{tool} {flags:?}");
            }
        }
    }
}

/// A stand-in sorter that records its arguments and the name it was found under.
fn recording_sorter(dir: &Path, tag: &str) {
    use std::os::unix::fs::PermissionsExt;
    let path = dir.join("tag2sort");
    std::fs::write(
        &path,
        format!("#!/bin/sh\necho \"{tag} $*\" >> \"$SORTER_LOG\"\n"),
    )
    .unwrap();
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
}

/// Runs a copied `SNPsplit` by bare name, with a decoy `tag2sort` earlier on `PATH`.
fn run_with_sorters(extra: &[&str], env_parallel: Option<&str>) -> String {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test/fixtures/se_basic");
    let work = TempDir::new().unwrap();
    let bin = work.path().join("bin");
    let decoy = work.path().join("decoy");
    std::fs::create_dir_all(&bin).unwrap();
    std::fs::create_dir_all(&decoy).unwrap();
    std::fs::copy(cargo_bin("snpsplit"), bin.join("SNPsplit")).unwrap();
    recording_sorter(&bin, "own");
    recording_sorter(&decoy, "decoy");
    std::fs::copy(fixture.join("input.sam"), work.path().join("input.sam")).unwrap();
    let log = work.path().join("sorter.log");
    let path = format!(
        "{}:{}:{}",
        decoy.display(),
        bin.display(),
        std::env::var("PATH").unwrap()
    );

    let mut cmd = Command::new("SNPsplit");
    cmd.current_dir(work.path())
        .env("PATH", path)
        .env("SORTER_LOG", &log)
        .env("SNPSPLIT_NO_SLEEP", "1")
        .env_remove("SNPSPLIT_PARALLEL")
        .args([
            "--snp_file",
            fixture.join("snps.txt").to_str().unwrap(),
            "--single_end",
            "--sam",
        ])
        .args(extra)
        .arg("input.sam");
    if let Some(n) = env_parallel {
        cmd.env("SNPSPLIT_PARALLEL", n);
    }
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::read_to_string(log).unwrap()
}

#[test]
fn the_sorter_is_the_one_next_to_the_running_binary_and_gets_the_worker_count() {
    let calls = run_with_sorters(&["--parallel", "3"], None);
    assert!(calls.starts_with("own "), "{calls}");
    assert!(calls.contains("--parallel 3"), "{calls}");
}

#[test]
fn the_worker_count_from_the_environment_reaches_the_sorter_too() {
    let calls = run_with_sorters(&[], Some("5"));
    assert!(calls.contains("--parallel 5"), "{calls}");
}

#[test]
fn a_symlinked_snpsplit_runs_the_tag2sort_beside_its_target() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test/fixtures/se_basic");
    let work = TempDir::new().unwrap();
    let (bin, elsewhere) = (work.path().join("bin"), work.path().join("elsewhere"));
    std::fs::create_dir_all(&bin).unwrap();
    std::fs::create_dir_all(&elsewhere).unwrap();
    std::fs::copy(cargo_bin("snpsplit"), bin.join("SNPsplit")).unwrap();
    recording_sorter(&bin, "own");
    symlink(bin.join("SNPsplit"), elsewhere.join("SNPsplit")).unwrap();
    std::fs::copy(fixture.join("input.sam"), work.path().join("input.sam")).unwrap();
    let log = work.path().join("sorter.log");

    let out = Command::new(elsewhere.join("SNPsplit"))
        .current_dir(work.path())
        .env("SORTER_LOG", &log)
        .env("SNPSPLIT_NO_SLEEP", "1")
        .args([
            "--snp_file",
            fixture.join("snps.txt").to_str().unwrap(),
            "--single_end",
            "--sam",
            "input.sam",
        ])
        .output()
        .unwrap();

    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(std::fs::read_to_string(log).unwrap().starts_with("own "));
}

#[test]
fn the_subcommand_form_sorts_without_a_tag2sort_beside_it() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test/fixtures/se_basic");
    let work = TempDir::new().unwrap();
    let alone = work.path().join("snpsplit");
    std::fs::copy(cargo_bin("snpsplit"), &alone).unwrap();
    std::fs::copy(fixture.join("input.sam"), work.path().join("input.sam")).unwrap();

    let out = Command::new(&alone)
        .current_dir(work.path())
        .env("SNPSPLIT_NO_SLEEP", "1")
        .args([
            "tag",
            "--snp_file",
            fixture.join("snps.txt").to_str().unwrap(),
            "--single_end",
            "--sam",
            "input.sam",
        ])
        .output()
        .unwrap();

    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        work.path().join("input.genome1.sam").exists()
            || work.path().join("input.genome1.bam").exists()
    );
}
