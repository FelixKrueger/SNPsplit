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

    assert!(out.status.success(), "exit status was {:?}", out.status.code());
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
