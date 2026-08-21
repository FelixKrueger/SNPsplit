//! Fetching the reference inputs the genome preparation needs.
//!
//! Running `SNPsplit_genome_preparation` means first finding and downloading two things by
//! hand: the Mouse Genomes Project SNP VCF and a folder of per-chromosome FastA files. Both
//! URLs live in comments in the Perl source and in the documentation, which is where users go
//! looking for them. A tool that already knows which build it wants can fetch them.
//!
//! This is the one thing in the port that the Perl does not do, so it is held to a different
//! standard: nothing here runs unless `--download` was given, and nothing here overwrites a
//! path the user named.

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use sha2::{Digest, Sha256};

/// Where the two inputs come from.
pub struct Sources {
    pub vcf_url: String,
    pub vcf_name: String,
    pub genome_base: String,
}

/// The MGP release to fetch, and the Ensembl release the genome comes from.
///
/// v8 is the current release and the default; v7 is the older combined SNP and INDEL file,
/// selected by `--v7_VCF`, and lives on a different host under a different name.
pub fn sources(v7: bool, genome_build: &str, ensembl_release: u32) -> Sources {
    if v7 {
        Sources {
            vcf_url: "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz".to_string(),
            vcf_name: "mgp_REL2005_snps_indels.vcf.gz".to_string(),
            genome_base: format!(
                "https://ftp.ensembl.org/pub/release-{ensembl_release}/fasta/mus_musculus/dna"
            ),
        }
    } else {
        Sources {
            vcf_url: "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz".to_string(),
            vcf_name: "mgp_REL2021_snps.vcf.gz".to_string(),
            genome_base: format!(
                "https://ftp.ensembl.org/pub/release-{ensembl_release}/fasta/mus_musculus/dna"
            ),
        }
    }
    .with_build(genome_build)
}

impl Sources {
    fn with_build(self, _genome_build: &str) -> Self {
        // The build is carried by the Ensembl release rather than by the path, so nothing to
        // substitute; kept as a seam for when a second build is supported.
        self
    }
}

/// The chromosomes a mouse genome folder needs.
pub fn mouse_chromosomes() -> Vec<String> {
    let mut names: Vec<String> = (1..=19).map(|n| n.to_string()).collect();
    names.extend(["X", "Y", "MT"].iter().map(|s| s.to_string()));
    names
}

/// Fetch `url` into `dest`, resuming a partial file rather than restarting it.
///
/// A partial download is continued with a range request. The VCF is large enough that losing
/// one to a dropped connection is a real cost, and a half-file left behind by a previous run
/// is exactly the case worth handling.
pub fn fetch(url: &str, dest: &Path, err: &mut impl Write) -> Result<u64> {
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let existing = std::fs::metadata(dest).map(|m| m.len()).unwrap_or(0);

    let mut request = ureq::get(url);
    if existing > 0 {
        writeln!(err, "Resuming '{}' from {existing} bytes", dest.display())?;
        request = request.header("Range", &format!("bytes={existing}-"));
    } else {
        writeln!(err, "Fetching {url}")?;
    }

    let mut response = request
        .call()
        .with_context(|| format!("Failed to fetch {url}"))?;

    let status = response.status().as_u16();
    // 206 means the server honoured the range; 200 means it ignored it and is sending the
    // whole file, in which case what is on disk has to be discarded rather than appended to.
    let append = status == 206;
    if !(status == 200 || status == 206) {
        anyhow::bail!("Failed to fetch {url}: HTTP {status}\n");
    }

    let mut file = if append {
        let mut f = std::fs::OpenOptions::new().append(true).open(dest)?;
        f.seek(SeekFrom::End(0))?;
        f
    } else {
        File::create(dest).with_context(|| format!("Failed to write to {}", dest.display()))?
    };

    let mut reader = response.body_mut().as_reader();
    let written = std::io::copy(&mut reader, &mut file)
        .with_context(|| format!("Failed while writing {}", dest.display()))?;
    file.flush()?;

    Ok(if append { existing + written } else { written })
}

/// The SHA-256 of a file, as lowercase hex.
pub fn sha256(path: &Path) -> Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    // sha2 0.11 returns a byte array rather than something with a hex Display.
    Ok(hasher
        .finalize()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect())
}

/// Check that a downloaded VCF is complete and is a VCF.
///
/// The MGP release publishes no checksum, so the file is verified by decompressing the whole
/// gzip stream and requiring a VCF header. That catches the two failures that actually
/// happen: a truncated transfer, and an HTML error page saved under a `.vcf.gz` name.
pub fn verify_vcf(path: &Path) -> Result<()> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(flate2::read::MultiGzDecoder::new(file));

    let mut first = String::new();
    reader
        .read_line(&mut first)
        .with_context(|| format!("'{}' is not readable as gzip", path.display()))?;
    if !first.starts_with("##fileformat=VCF") {
        anyhow::bail!(
            "'{}' does not begin with a VCF header; the download is not a VCF\n",
            path.display()
        );
    }

    // Read to the end so a truncation anywhere in the file surfaces here rather than during
    // the run that uses it.
    let mut sink = std::io::sink();
    std::io::copy(&mut reader, &mut sink)
        .with_context(|| format!("'{}' is truncated or corrupt", path.display()))?;
    Ok(())
}

/// One entry of the record a download leaves behind.
pub struct Fetched {
    pub url: String,
    pub path: PathBuf,
    pub bytes: u64,
    pub sha256: String,
}

/// Append what was fetched to `<download_dir>/manifest.txt`.
///
/// So a later run can say whether the remote file has changed rather than silently using a
/// different input, which for a SNP annotation is the difference between two experiments.
pub fn record(dir: &Path, entries: &[Fetched]) -> Result<()> {
    let path = dir.join("manifest.txt");
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .with_context(|| format!("Failed to write to {}", path.display()))?;
    for entry in entries {
        writeln!(
            file,
            "{}\t{}\t{}\t{}",
            entry.path.display(),
            entry.bytes,
            entry.sha256,
            entry.url
        )?;
    }
    Ok(())
}

/// Parse an Ensembl `CHECKSUMS` file into name and checksum-line pairs.
///
/// Ensembl publishes BSD `sum` output: two numbers then the file name.
pub fn parse_checksums(text: &str) -> Vec<(String, String)> {
    text.lines()
        .filter_map(|line| {
            let mut parts = line.split_whitespace();
            let a = parts.next()?;
            let b = parts.next()?;
            let name = parts.next()?;
            Some((name.to_string(), format!("{a} {b}")))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_v8_release_is_the_default_and_v7_lives_elsewhere() {
        let v8 = sources(false, "GRCm39", 115);
        assert!(v8.vcf_url.contains("REL-2112-v8"));
        assert_eq!(v8.vcf_name, "mgp_REL2021_snps.vcf.gz");

        let v7 = sources(true, "GRCm39", 115);
        assert!(v7.vcf_url.contains("REL-2004-v7"));
        assert_eq!(v7.vcf_name, "mgp_REL2005_snps_indels.vcf.gz");
    }

    #[test]
    fn the_ensembl_release_appears_in_the_genome_url() {
        assert!(
            sources(false, "GRCm39", 115)
                .genome_base
                .contains("release-115")
        );
    }

    #[test]
    fn a_mouse_genome_has_twenty_two_sequences() {
        let names = mouse_chromosomes();
        assert_eq!(names.len(), 22);
        assert_eq!(names[0], "1");
        assert_eq!(names[19], "X");
        assert_eq!(names[21], "MT");
    }

    #[test]
    fn checksums_are_name_and_sum_pairs() {
        let parsed = parse_checksums(
            "12345 6789 Mus_musculus.GRCm39.dna.chromosome.1.fa.gz\n1 2 other.fa.gz\n",
        );
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].0, "Mus_musculus.GRCm39.dna.chromosome.1.fa.gz");
        assert_eq!(parsed[0].1, "12345 6789");
    }

    #[test]
    fn a_download_that_is_not_a_vcf_is_rejected() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("not.vcf.gz");
        let file = File::create(&path).unwrap();
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(b"<html>404 Not Found</html>\n").unwrap();
        gz.finish().unwrap();

        let Err(e) = verify_vcf(&path) else {
            panic!("an HTML error page saved as a VCF should not verify");
        };
        assert!(e.to_string().contains("does not begin with a VCF header"));
    }

    #[test]
    fn a_complete_vcf_verifies() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("ok.vcf.gz");
        let file = File::create(&path).unwrap();
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(b"##fileformat=VCFv4.2\n#CHROM\tPOS\n")
            .unwrap();
        gz.finish().unwrap();

        verify_vcf(&path).unwrap();
    }

    #[test]
    fn hashing_is_stable_and_content_dependent() {
        let dir = tempfile::TempDir::new().unwrap();
        let a = dir.path().join("a");
        File::create(&a).unwrap().write_all(b"abc").unwrap();
        // The SHA-256 of "abc", which is a published test vector rather than whatever this
        // implementation happens to produce.
        assert_eq!(
            sha256(&a).unwrap(),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }
}
