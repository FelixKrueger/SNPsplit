//! Helpers shared by the integration tests.

/// Whether samtools runs. Panics instead of returning false when SNPSPLIT_REQUIRE_SAMTOOLS is 1 or true.
pub fn samtools_present() -> bool {
    let present = std::process::Command::new("samtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success());
    let required = std::env::var("SNPSPLIT_REQUIRE_SAMTOOLS");
    if !present && matches!(required.as_deref(), Ok("1" | "true")) {
        panic!("samtools not found and SNPSPLIT_REQUIRE_SAMTOOLS is set");
    }
    present
}
