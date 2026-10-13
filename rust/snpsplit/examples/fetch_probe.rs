//! Fetches one small real file, to prove the HTTPS path works end to end.
//! Not part of the shipped binary; run by hand.
fn main() -> anyhow::Result<()> {
    let url = std::env::args().nth(1).expect("url");
    let dir = tempfile::tempdir()?;
    let dest = dir.path().join("probe");
    let bytes = snpsplit::genome_prep::download::fetch(&url, &dest, &mut std::io::stderr())?;
    println!(
        "fetched {bytes} bytes, sha256 {}",
        snpsplit::genome_prep::download::sha256(&dest)?
    );
    Ok(())
}
