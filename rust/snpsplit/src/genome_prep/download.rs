//! Fetching the reference inputs the genome preparation needs.
//!
//! Running `SNPsplit_genome_preparation` means first finding and downloading two things by
//! hand: the Mouse Genomes Project SNP VCF and a folder of per-chromosome FastA files. A tool
//! that already knows which build it wants can fetch them.
//!
//! Nothing here runs unless `--download` was given. A file or folder the user named is changed
//! only when it is recognisably the output of `--download`. Every transfer goes to
//! `<target>.part` and is renamed only once it has been checked, so a run started again after an
//! interruption finishes what is missing.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use md5::Md5;
use sha2::{Digest, Sha256};

/// File name of the MGP v8 SNP VCF, the only VCF `--download` fetches.
pub const V8_VCF_NAME: &str = "mgp_REL2021_snps.vcf.gz";

/// Written into a genome folder that `--download` fills; holds the build and Ensembl release.
pub const MARKER: &str = ".snpsplit_download";

/// Default Ensembl release for each build `--download` can fetch.
const GRCM39_RELEASE: u32 = 115;
const GRCM38_RELEASE: u32 = 102;

/// Where the two inputs come from.
pub struct Sources {
    pub vcf_url: String,
    pub vcf_md5_url: String,
    pub genome_base: String,
}

/// The MGP v8 VCF (GRCm39) and the Ensembl genome folder of `ensembl_release`.
pub fn sources(ensembl_release: u32) -> Sources {
    let base = "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels";
    Sources {
        vcf_url: format!("{base}/{V8_VCF_NAME}"),
        vcf_md5_url: format!("{base}/{V8_VCF_NAME}.md5"),
        genome_base: format!(
            "https://ftp.ensembl.org/pub/release-{ensembl_release}/fasta/mus_musculus/dna"
        ),
    }
}

/// The chromosomes a mouse genome folder needs.
pub fn mouse_chromosomes() -> Vec<String> {
    let mut names: Vec<String> = (1..=19).map(|n| n.to_string()).collect();
    names.extend(["X", "Y", "MT"].iter().map(|s| s.to_string()));
    names
}

/// Whether a path is fetched by `--download` or used as the user gave it.
#[derive(Debug, PartialEq, Eq)]
pub enum Owner {
    Ours,
    Theirs,
}

/// A VCF path is ours unless the user named a file that exists and that we did not fetch.
pub fn vcf_owner(path: &Path, chosen_by_download: bool, manifest: &Manifest) -> Owner {
    let v8_name = path.file_name().is_some_and(|n| n == V8_VCF_NAME);
    if chosen_by_download || !path.exists() || manifest.contains(path) || v8_name {
        Owner::Ours
    } else {
        Owner::Theirs
    }
}

/// A genome folder is ours unless the user named one that holds FastA files we did not fetch.
pub fn genome_owner(dir: &Path, chosen_by_download: bool, build: &str) -> Owner {
    if chosen_by_download
        || !dir.is_dir()
        || !has_fasta(dir)
        || dir.join(MARKER).exists()
        || left_by_1_0_0(dir, build)
    {
        Owner::Ours
    } else {
        Owner::Theirs
    }
}

/// Whether a folder holds at least one FastA file.
pub fn has_fasta(dir: &Path) -> bool {
    names_in(dir).any(|n| n.ends_with(".fa") || n.ends_with(".fasta"))
}

/// Whether a folder looks like an interrupted 1.0.0 download of `build`: a chromosome archive of
/// that build, and FastA files named only `<chr>.fa`.
fn left_by_1_0_0(dir: &Path, build: &str) -> bool {
    let archive = format!("Mus_musculus.{build}.dna.chromosome.");
    let chromosomes: Vec<String> = mouse_chromosomes()
        .iter()
        .map(|c| format!("{c}.fa"))
        .collect();
    let names: Vec<String> = names_in(dir).collect();
    names
        .iter()
        .any(|n| n.starts_with(&archive) && (n.ends_with(".fa.gz") || n.ends_with(".fa.gz.part")))
        && names
            .iter()
            .filter(|n| n.ends_with(".fa") || n.ends_with(".fasta"))
            .all(|n| chromosomes.contains(n))
}

/// The Ensembl release a folder was filled from, if `--download` filled it with `build`.
pub fn marked_release(dir: &Path, build: &str) -> Option<u32> {
    let text = std::fs::read_to_string(dir.join(MARKER)).ok()?;
    let (marked, release) = text.trim().split_once(' ')?;
    (marked == build).then(|| release.parse().ok()).flatten()
}

fn names_in(dir: &Path) -> impl Iterator<Item = String> {
    std::fs::read_dir(dir)
        .into_iter()
        .flatten()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
}

/// Refuse a download that cannot work, before anything is fetched; returns the Ensembl release.
///
/// `named_vcf` is the `--vcf_file` the user gave. The VCF is checked against the genome build
/// only when both come from MGP and Ensembl; a user's own genome is their call.
pub fn check_request(
    v7: bool,
    named_vcf: Option<&Path>,
    fetch_vcf: bool,
    fetch_genome: bool,
    build: &str,
    release: Option<u32>,
) -> Result<u32> {
    if v7 && fetch_vcf {
        match named_vcf {
            Some(p) if p.exists() => anyhow::bail!(
                "--v7_VCF was given, but '{}' is the v8 VCF. Remove --v7_VCF, or give the v7 file with --vcf_file.\n",
                p.display()
            ),
            Some(p) => anyhow::bail!(
                "'{}' does not exist, and the v7 VCF can no longer be downloaded. Check the path given with --vcf_file.\n",
                p.display()
            ),
            None => anyhow::bail!(
                "The v7 VCF (mgp_REL2005_snps_indels.vcf.gz) can no longer be downloaded. Supply it with --vcf_file.\n"
            ),
        }
    }
    if !fetch_genome {
        return Ok(release.unwrap_or(GRCM39_RELEASE));
    }

    let (default, valid, range) = match build {
        "GRCm39" => (
            GRCM39_RELEASE,
            release.is_none_or(|r| r >= 103),
            "103 or later",
        ),
        "GRCm38" => (
            GRCM38_RELEASE,
            release.is_none_or(|r| (76..=102).contains(&r)),
            "76 to 102",
        ),
        other => anyhow::bail!(
            "--download can fetch the GRCm39 or GRCm38 genome, not '{other}'. Set --genome_build, or supply the genome with --reference_genome.\n"
        ),
    };
    if !valid {
        anyhow::bail!(
            "Ensembl has {build} in releases {range}, not release {}. Set --ensembl_release accordingly.\n",
            release.unwrap_or(default)
        );
    }

    let v8_vcf =
        fetch_vcf || named_vcf.is_some_and(|p| p.file_name().is_some_and(|n| n == V8_VCF_NAME));
    if v8_vcf && build != "GRCm39" {
        anyhow::bail!(
            "The MGP v8 VCF uses GRCm39 coordinates, but the genome to download is {build}. Use --genome_build GRCm39.\n"
        );
    }
    if v7 && build != "GRCm38" {
        anyhow::bail!(
            "The v7 VCF uses GRCm38 coordinates, but the genome to download is {build}. Use --genome_build GRCm38 (Ensembl release 102 is used).\n"
        );
    }
    Ok(release.unwrap_or(default))
}

/// Fetch `url` into `dest`, resuming a partial file; returns the size of `dest` afterwards.
///
/// A 416 for a resumed file means the range starts at or past the end of the remote file. That
/// is taken as "nothing left to transfer"; the caller's check decides whether the file is right.
pub fn fetch(url: &str, dest: &Path, err: &mut impl Write) -> Result<u64> {
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create {}", parent.display()))?;
    }
    let existing = std::fs::metadata(dest).map(|m| m.len()).unwrap_or(0);

    let mut request = agent().get(url);
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
    let append = match status {
        206 if existing > 0 => true,
        200 => {
            if existing > 0 {
                writeln!(
                    err,
                    "The server sent the whole file; starting '{}' again",
                    dest.display()
                )?;
            }
            false
        }
        416 if existing > 0 => return Ok(existing),
        _ => anyhow::bail!("Failed to fetch {url}: HTTP {status}\n"),
    };
    let mut file = if append {
        std::fs::OpenOptions::new().append(true).open(dest)
    } else {
        File::create(dest)
    }
    .with_context(|| format!("Failed to write to {}", dest.display()))?;

    let mut reader = response.body_mut().as_reader();
    let written = std::io::copy(&mut reader, &mut file).with_context(|| {
        format!(
            "Failed while fetching {url} into {}; run the same command again to resume",
            dest.display()
        )
    })?;
    file.flush()
        .with_context(|| format!("Failed to write to {}", dest.display()))?;
    Ok(if append { existing + written } else { written })
}

/// An agent that reports 4xx and 5xx as responses, so `fetch` can tell a 416 from a 404.
fn agent() -> ureq::Agent {
    ureq::Agent::new_with_config(
        ureq::Agent::config_builder()
            .http_status_as_error(false)
            .build(),
    )
}

/// The body of a small text resource, such as a published checksum.
fn get_text(url: &str) -> Result<String> {
    let mut response = agent()
        .get(url)
        .call()
        .with_context(|| format!("Failed to fetch {url}"))?;
    let status = response.status().as_u16();
    if status != 200 {
        anyhow::bail!("Failed to fetch {url}: HTTP {status}\n");
    }
    response
        .body_mut()
        .read_to_string()
        .with_context(|| format!("Failed to read {url}"))
}

/// The MD5 of `data`, as lowercase hex.
pub fn md5_hex(data: &[u8]) -> String {
    hex(&Md5::digest(data))
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// MD5 and SHA-256 of a file in one read, as lowercase hex.
fn hashes(path: &Path) -> Result<(String, String)> {
    let mut file =
        File::open(path).with_context(|| format!("Failed to read {}", path.display()))?;
    let mut md5 = Md5::new();
    let mut sha = Sha256::new();
    let mut buffer = vec![0u8; 1 << 20];
    loop {
        let read = file
            .read(&mut buffer)
            .with_context(|| format!("Failed to read {}", path.display()))?;
        if read == 0 {
            break;
        }
        md5.update(&buffer[..read]);
        sha.update(&buffer[..read]);
    }
    Ok((hex(&md5.finalize()), hex(&sha.finalize())))
}

/// The SHA-256 of a file, as lowercase hex.
pub fn sha256(path: &Path) -> Result<String> {
    Ok(hashes(path)?.1)
}

/// Whether a file starts with the gzip signature.
fn is_gzip(path: &Path) -> bool {
    let mut magic = [0u8; 2];
    File::open(path)
        .and_then(|mut f| f.read_exact(&mut magic))
        .is_ok()
        && magic == [0x1f, 0x8b]
}

/// Whether the first line of a gzipped file is a VCF header.
pub fn has_vcf_header(path: &Path) -> bool {
    let Ok(file) = File::open(path) else {
        return false;
    };
    let mut first = String::new();
    BufReader::new(flate2::read::MultiGzDecoder::new(file))
        .read_line(&mut first)
        .is_ok()
        && first.starts_with("##fileformat=VCF")
}

/// Check a VCF against the published MD5; returns its SHA-256 for the manifest.
fn check_vcf(path: &Path, expected_md5: &str) -> Result<String> {
    if !has_vcf_header(path) {
        anyhow::bail!("'{}' does not begin with a VCF header\n", path.display());
    }
    let (md5, sha) = hashes(path)?;
    if md5 != expected_md5 {
        anyhow::bail!(
            "'{}' does not match the published MD5 ({md5}, expected {expected_md5})\n",
            path.display()
        );
    }
    Ok(sha)
}

fn published_md5(src: &Sources) -> Result<String> {
    let text = get_text(&src.vcf_md5_url)?;
    text.split_whitespace()
        .next()
        .filter(|h| h.len() == 32 && h.chars().all(|c| c.is_ascii_hexdigit()))
        .map(|h| h.to_ascii_lowercase())
        .with_context(|| format!("{} does not hold an MD5\n", src.vcf_md5_url))
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

fn rename(from: &Path, to: &Path) -> Result<()> {
    std::fs::rename(from, to)
        .with_context(|| format!("Failed to rename {} to {}", from.display(), to.display()))
}

fn remove(path: &Path) -> Result<()> {
    std::fs::remove_file(path).with_context(|| format!("Failed to remove {}", path.display()))
}

fn remove_if_present(path: &Path) -> Result<()> {
    if path.exists() { remove(path) } else { Ok(()) }
}

/// Fetch the v8 VCF into `target`, finishing whatever an earlier run left.
///
/// `named` is true when the user gave `target` with `--vcf_file`. An existing named file that
/// fails the check is left as it is, unless an earlier `--download` wrote it.
pub fn fetch_vcf(
    target: &Path,
    named: bool,
    src: &Sources,
    manifest: &mut Manifest,
    err: &mut impl Write,
) -> Result<()> {
    let part = with_suffix(target, ".part");

    if target.exists() && manifest.lists(target) {
        remove_if_present(&part)?;
        writeln!(err, "Using '{}' (listed in manifest.txt)", target.display())?;
        return Ok(());
    }
    let expected = published_md5(src)?;
    let unusable = || format!("The VCF fetched from {} is not usable", src.vcf_url);

    if target.exists() {
        if let Ok(sha) = check_vcf(target, &expected) {
            remove_if_present(&part)?;
            writeln!(err, "Using '{}' (MD5 matches)", target.display())?;
            return manifest.add(target, &src.vcf_url, &sha);
        }
        let ours = !target.is_symlink() && (!named || manifest.contains(target));
        if !ours {
            anyhow::bail!(
                "'{}' is not the MGP v8 VCF published at {}. It was left unchanged: remove it to download it again, or give another file with --vcf_file.\n",
                target.display(),
                src.vcf_url
            );
        }
        if is_gzip(target) {
            rename(target, &part)?;
        } else {
            remove(target)?;
        }
    }

    let mut resumed = part.exists();
    if resumed && !is_gzip(&part) {
        writeln!(err, "Restarting '{}' (not a gzip file)", part.display())?;
        remove(&part)?;
        resumed = false;
    }

    fetch(&src.vcf_url, &part, err)?;
    let sha = match check_vcf(&part, &expected) {
        Ok(sha) => sha,
        Err(_) if resumed => {
            writeln!(
                err,
                "Restarting '{}' (the resumed file did not check out)",
                part.display()
            )?;
            remove(&part)?;
            fetch(&src.vcf_url, &part, err)?;
            check_vcf(&part, &expected)
                .inspect_err(|_| {
                    let _ = std::fs::remove_file(&part);
                })
                .with_context(unusable)?
        }
        Err(e) => {
            let _ = std::fs::remove_file(&part);
            return Err(e.context(unusable()));
        }
    };

    rename(&part, target)?;
    writeln!(err, "Verified '{}'", target.display())?;
    manifest.add(target, &src.vcf_url, &sha)
}

/// Fetch every chromosome of `build` that `dir` does not yet hold.
///
/// A fresh archive that does not unpack stops the run; one that was resumed is fetched once more
/// from scratch first. Either way at most one resume and one restart happen per chromosome.
pub fn fetch_genome(
    dir: &Path,
    base_url: &str,
    build: &str,
    release: u32,
    chromosomes: &[String],
    manifest: &mut Manifest,
    err: &mut impl Write,
) -> Result<()> {
    let archive = |chr: &str| dir.join(format!("Mus_musculus.{build}.dna.chromosome.{chr}.fa.gz"));
    // A chromosome is finished once its archive is gone: it is removed only after unpacking.
    let missing: Vec<&String> = chromosomes
        .iter()
        .filter(|chr| {
            let gz = archive(chr);
            !dir.join(format!("{chr}.fa")).exists()
                || gz.exists()
                || with_suffix(&gz, ".part").exists()
        })
        .collect();
    if missing.is_empty() {
        writeln!(
            err,
            "Using the reference genome already present at '{}'",
            dir.display()
        )?;
        return Ok(());
    }
    std::fs::create_dir_all(dir).with_context(|| format!("Failed to create {}", dir.display()))?;
    claim(dir, build, release)?;

    for chr in missing {
        let gz = archive(chr);
        let name = gz
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let url = format!("{base_url}/{name}");
        let target = dir.join(format!("{chr}.fa"));
        let gz_part = with_suffix(&gz, ".part");

        if gz.exists() {
            rename(&gz, &gz_part)?;
        }
        let mut resumed = gz_part.exists();
        if resumed && !is_gzip(&gz_part) {
            writeln!(err, "Restarting '{}' (not a gzip file)", gz_part.display())?;
            remove(&gz_part)?;
            resumed = false;
        }

        let mut bytes = fetch(&url, &gz_part, err)?;
        rename(&gz_part, &gz)?;
        match unpack(&gz, &target) {
            Ok(()) => {}
            Err(Unpack::Corrupt(e)) if resumed => {
                writeln!(err, "Restarting '{}' ({e:#})", gz.display())?;
                remove(&gz)?;
                bytes = fetch(&url, &gz_part, err)?;
                rename(&gz_part, &gz)?;
                match unpack(&gz, &target) {
                    Ok(()) => {}
                    Err(Unpack::Corrupt(e)) => {
                        return Err(
                            e.context(format!("The archive fetched from {url} is not usable"))
                        );
                    }
                    Err(Unpack::Write(e)) => return Err(e),
                }
            }
            Err(Unpack::Corrupt(e)) => {
                return Err(e.context(format!("The archive fetched from {url} is not usable")));
            }
            Err(Unpack::Write(e)) => return Err(e),
        }

        let (_, sha) = hashes(&gz)?;
        manifest.add_entry(&target, bytes, &sha, &url)?;
        remove(&gz)?;
    }
    Ok(())
}

/// Mark `dir` as filled by `--download` from `build` at `release`, or refuse if it holds another.
fn claim(dir: &Path, build: &str, release: u32) -> Result<()> {
    let marker = dir.join(MARKER);
    let wanted = format!("{build} {release}");
    let found = std::fs::read_to_string(&marker).unwrap_or_default();
    let found = found.trim();
    if !found.is_empty() && found != wanted {
        anyhow::bail!(
            "'{}' was filled from Ensembl {found}, not {wanted}. Use the same --genome_build and --ensembl_release, or another folder.\n",
            dir.display()
        );
    }
    std::fs::write(&marker, format!("{wanted}\n"))
        .with_context(|| format!("Failed to write {}", marker.display()))
}

/// Why a chromosome archive did not unpack: its content, or reading or writing a file.
enum Unpack {
    Corrupt(anyhow::Error),
    Write(anyhow::Error),
}

/// Unpack a chromosome archive to `<target>.part`, then rename it to `target`.
fn unpack(source: &Path, target: &Path) -> std::result::Result<(), Unpack> {
    let part = with_suffix(target, ".part");
    let write_failed = |e: std::io::Error| {
        Unpack::Write(anyhow::Error::new(e).context(format!("Failed to write {}", part.display())))
    };
    let input = File::open(source).map_err(|e| {
        Unpack::Write(anyhow::Error::new(e).context(format!("Failed to read {}", source.display())))
    })?;
    let mut reader = flate2::read::MultiGzDecoder::new(input);
    let mut output = File::create(&part).map_err(write_failed)?;
    let mut buffer = vec![0u8; 1 << 20];
    loop {
        let read = match reader.read(&mut buffer) {
            Ok(0) => break,
            Ok(n) => n,
            Err(e) => {
                let _ = std::fs::remove_file(&part);
                return Err(Unpack::Corrupt(anyhow::Error::new(e).context(format!(
                    "'{}' is truncated or corrupt",
                    source.display()
                ))));
            }
        };
        if let Err(e) = output.write_all(&buffer[..read]) {
            let _ = std::fs::remove_file(&part);
            return Err(write_failed(e));
        }
    }
    if let Err(e) = output.flush() {
        let _ = std::fs::remove_file(&part);
        return Err(write_failed(e));
    }
    rename(&part, target).map_err(Unpack::Write)
}

/// What `--download` has fetched into a folder: `manifest.txt`, one tab-separated line per file
/// (path, bytes, sha256, url). Bytes and SHA-256 are of what was downloaded, so for a chromosome
/// they describe the `.fa.gz`, not the unpacked `.fa`. Later lines for the same path win.
pub struct Manifest {
    file: PathBuf,
    entries: HashMap<PathBuf, u64>,
}

impl Manifest {
    /// Read `<dir>/manifest.txt`; a missing file is an empty manifest.
    pub fn load(dir: &Path) -> Result<Self> {
        let file = dir.join("manifest.txt");
        let mut entries = HashMap::new();
        let text = match std::fs::read_to_string(&file) {
            Ok(text) => text,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => String::new(),
            Err(e) => {
                return Err(
                    anyhow::Error::new(e).context(format!("Failed to read {}", file.display()))
                );
            }
        };
        for line in text.lines() {
            let mut fields = line.split('\t');
            let (Some(path), Some(bytes)) = (fields.next(), fields.next()) else {
                continue;
            };
            if let Ok(bytes) = bytes.parse() {
                entries.insert(key(Path::new(path)), bytes);
            }
        }
        Ok(Self { file, entries })
    }

    /// Whether the manifest has any entry for `path`.
    pub fn contains(&self, path: &Path) -> bool {
        self.entries.contains_key(&key(path))
    }

    /// Whether `path` is listed and its size on disk still matches; for files stored as fetched.
    pub fn lists(&self, path: &Path) -> bool {
        let size = std::fs::metadata(path).map(|m| m.len()).ok();
        self.entries
            .get(&key(path))
            .is_some_and(|&b| Some(b) == size)
    }

    fn add(&mut self, path: &Path, url: &str, sha256: &str) -> Result<()> {
        let bytes = std::fs::metadata(path)
            .with_context(|| format!("Failed to read {}", path.display()))?
            .len();
        self.add_entry(path, bytes, sha256, url)
    }

    fn add_entry(&mut self, path: &Path, bytes: u64, sha256: &str, url: &str) -> Result<()> {
        if let Some(parent) = self.file.parent() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("Failed to create {}", parent.display()))?;
        }
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.file)
            .with_context(|| format!("Failed to write to {}", self.file.display()))?;
        let path = key(path);
        writeln!(file, "{}\t{bytes}\t{sha256}\t{url}", path.display())
            .with_context(|| format!("Failed to write to {}", self.file.display()))?;
        self.entries.insert(path, bytes);
        Ok(())
    }
}

/// A path in the form manifest entries are compared in: absolute, symlinks resolved.
fn key(path: &Path) -> PathBuf {
    if let Ok(p) = path.canonicalize() {
        return p;
    }
    let absolute = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    match (absolute.parent(), absolute.file_name()) {
        (Some(parent), Some(name)) => parent
            .canonicalize()
            .map(|p| p.join(name))
            .unwrap_or(absolute.clone()),
        _ => absolute,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_gz(path: &Path, data: &[u8]) {
        let mut gz = flate2::write::GzEncoder::new(
            File::create(path).unwrap(),
            flate2::Compression::default(),
        );
        gz.write_all(data).unwrap();
        gz.finish().unwrap();
    }

    #[test]
    fn only_the_v8_vcf_is_fetched_and_its_md5_sits_beside_it() {
        let s = sources(115);
        assert!(
            s.vcf_url
                .ends_with("/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz")
        );
        assert_eq!(s.vcf_md5_url, format!("{}.md5", s.vcf_url));
        assert!(s.genome_base.contains("release-115"));
    }

    #[test]
    fn a_mouse_genome_has_twenty_two_sequences() {
        let names = mouse_chromosomes();
        assert_eq!(names.len(), 22);
        assert_eq!(names[19], "X");
        assert_eq!(names[21], "MT");
    }

    #[test]
    fn requests_that_cannot_work_are_refused_up_front() {
        let own = Path::new("own.vcf.gz");
        let v8 = Path::new("dir/mgp_REL2021_snps.vcf.gz");
        let v7 = Path::new("dir/mgp_REL2005_snps_indels.vcf.gz");
        assert!(check_request(true, None, true, false, "GRCm39", None).is_err());
        let named_v7 = check_request(true, Some(v7), true, false, "GRCm38", None).unwrap_err();
        assert!(
            named_v7.to_string().contains("does not exist"),
            "{named_v7}"
        );
        assert_eq!(
            check_request(true, Some(own), false, true, "GRCm38", None).unwrap(),
            102
        );
        assert!(check_request(true, Some(own), false, true, "GRCm39", None).is_err());
        assert_eq!(
            check_request(false, None, true, true, "GRCm39", None).unwrap(),
            115
        );
        assert_eq!(
            check_request(false, Some(own), false, true, "GRCm38", None).unwrap(),
            102
        );
        assert!(check_request(false, Some(own), false, true, "GRCm38", Some(115)).is_err());
        assert!(check_request(false, Some(own), false, true, "GRCm39", Some(102)).is_err());
        assert!(check_request(false, Some(own), false, true, "hg38", None).is_err());
        assert!(check_request(false, None, true, true, "GRCm38", None).is_err());
        assert!(check_request(false, Some(v8), false, true, "GRCm38", None).is_err());
        // The user's own genome is not checked against the VCF.
        assert_eq!(
            check_request(false, None, true, false, "GRCm38", Some(99)).unwrap(),
            99
        );
    }

    #[test]
    fn a_vcf_is_ours_unless_it_is_a_users_existing_file() {
        let dir = tempfile::TempDir::new().unwrap();
        let empty = Manifest::load(dir.path()).unwrap();
        let own = dir.path().join("own.vcf.gz");
        let named_v8 = dir.path().join(V8_VCF_NAME);
        std::fs::write(&own, b"x").unwrap();
        std::fs::write(&named_v8, b"x").unwrap();

        assert_eq!(
            vcf_owner(&dir.path().join("x.vcf.gz"), true, &empty),
            Owner::Ours
        );
        assert_eq!(
            vcf_owner(&dir.path().join("missing.vcf.gz"), false, &empty),
            Owner::Ours
        );
        assert_eq!(vcf_owner(&named_v8, false, &empty), Owner::Ours);
        assert_eq!(vcf_owner(&own, false, &empty), Owner::Theirs);

        let mut listed = Manifest::load(dir.path()).unwrap();
        listed.add(&own, "http://x", "0").unwrap();
        assert_eq!(vcf_owner(&own, false, &listed), Owner::Ours);
    }

    #[test]
    fn a_genome_folder_is_ours_unless_it_holds_a_users_fasta() {
        let dir = tempfile::TempDir::new().unwrap();
        let folder = |name: &str, files: &[&str]| {
            let d = dir.path().join(name);
            std::fs::create_dir(&d).unwrap();
            for f in files {
                std::fs::write(d.join(f), ">1\nA\n").unwrap();
            }
            d
        };
        let empty = folder("empty", &[]);
        let theirs = folder("theirs", &["1.fa"]);
        let marked = folder("marked", &["1.fa", MARKER]);
        let left = folder(
            "left",
            &["1.fa", "Mus_musculus.GRCm39.dna.chromosome.2.fa.gz"],
        );
        let left_part = folder(
            "left_part",
            &["1.fa", "Mus_musculus.GRCm39.dna.chromosome.2.fa.gz.part"],
        );
        let kept_archives = folder(
            "kept",
            &["chr1.fa", "Mus_musculus.GRCm39.dna.chromosome.1.fa.gz"],
        );
        let other_build = folder(
            "grcm38",
            &["1.fa", "Mus_musculus.GRCm38.dna.chromosome.2.fa.gz"],
        );
        let own = |d: &Path| genome_owner(d, false, "GRCm39");

        assert_eq!(genome_owner(&theirs, true, "GRCm39"), Owner::Ours);
        assert_eq!(own(&dir.path().join("missing")), Owner::Ours);
        assert_eq!(own(&empty), Owner::Ours);
        assert_eq!(own(&marked), Owner::Ours);
        assert_eq!(own(&left), Owner::Ours);
        assert_eq!(own(&left_part), Owner::Ours);
        assert_eq!(own(&theirs), Owner::Theirs);
        assert_eq!(
            own(&kept_archives),
            Owner::Theirs,
            "a user's FastA beside its archives"
        );
        assert_eq!(
            own(&other_build),
            Owner::Theirs,
            "an archive of another build"
        );
    }

    #[test]
    fn a_folder_from_another_release_is_not_mixed_into() {
        let dir = tempfile::TempDir::new().unwrap();
        claim(dir.path(), "GRCm39", 115).unwrap();
        claim(dir.path(), "GRCm39", 115).unwrap();
        let e = claim(dir.path(), "GRCm39", 114).unwrap_err();
        assert!(e.to_string().contains("GRCm39 115"), "{e}");
        assert_eq!(marked_release(dir.path(), "GRCm39"), Some(115));
        assert_eq!(marked_release(dir.path(), "GRCm38"), None);
    }

    #[test]
    fn the_v7_message_says_why_the_named_file_cannot_be_used() {
        let dir = tempfile::TempDir::new().unwrap();
        let v8 = dir.path().join(V8_VCF_NAME);
        std::fs::write(&v8, b"x").unwrap();
        let e = check_request(true, Some(&v8), true, false, "GRCm38", None).unwrap_err();
        assert!(e.to_string().contains("is the v8 VCF"), "{e}");
    }

    #[test]
    fn a_failed_write_while_unpacking_is_not_corrupt() {
        let dir = tempfile::TempDir::new().unwrap();
        let gz = dir.path().join("1.fa.gz");
        write_gz(&gz, b">1\nACGT\n");
        let target = dir.path().join("1.fa");
        std::fs::create_dir(with_suffix(&target, ".part")).unwrap();
        assert!(matches!(unpack(&gz, &target), Err(Unpack::Write(_))));
        assert!(!target.exists());
    }

    #[test]
    fn manifest_paths_match_however_they_are_spelled_and_the_last_entry_wins() {
        let dir = tempfile::TempDir::new().unwrap();
        let f = dir.path().join("a.fa");
        std::fs::write(&f, b"12345").unwrap();
        let line = |p: &Path, n: u64| format!("{}\t{n}\tsha\turl\n", p.display());
        let roundabout = dir.path().join(".").join("a.fa");
        std::fs::write(
            dir.path().join("manifest.txt"),
            line(&roundabout, 3) + &line(&f, 5),
        )
        .unwrap();

        let m = Manifest::load(dir.path()).unwrap();
        assert!(m.contains(&f) && m.contains(&roundabout));
        assert!(m.lists(&f), "the later, matching entry should win");

        std::fs::write(&f, b"1234567").unwrap();
        assert!(!m.lists(&f), "a size mismatch is not listed");
    }

    #[test]
    fn an_unreadable_manifest_is_an_error_not_an_empty_one() {
        let dir = tempfile::TempDir::new().unwrap();
        std::fs::create_dir(dir.path().join("manifest.txt")).unwrap();
        assert!(Manifest::load(dir.path()).is_err());
    }

    #[test]
    fn a_vcf_checks_out_only_with_a_header_and_the_published_md5() {
        let dir = tempfile::TempDir::new().unwrap();
        let ok = dir.path().join("ok.vcf.gz");
        write_gz(&ok, b"##fileformat=VCFv4.2\n#CHROM\tPOS\n");
        let md5 = md5_hex(&std::fs::read(&ok).unwrap());
        check_vcf(&ok, &md5).unwrap();
        assert!(check_vcf(&ok, &"0".repeat(32)).is_err());

        let page = dir.path().join("page.vcf.gz");
        write_gz(&page, b"<html>404 Not Found</html>\n");
        let md5 = md5_hex(&std::fs::read(&page).unwrap());
        assert!(
            check_vcf(&page, &md5)
                .unwrap_err()
                .to_string()
                .contains("VCF header")
        );
    }

    #[test]
    fn a_truncated_archive_is_corrupt_rather_than_a_write_failure() {
        let dir = tempfile::TempDir::new().unwrap();
        let gz = dir.path().join("1.fa.gz");
        write_gz(&gz, &b">1\nACGT\n".repeat(1000));
        let full = std::fs::read(&gz).unwrap();
        std::fs::write(&gz, &full[..full.len() / 2]).unwrap();

        let target = dir.path().join("1.fa");
        assert!(matches!(unpack(&gz, &target), Err(Unpack::Corrupt(_))));
        assert!(!target.exists() && !with_suffix(&target, ".part").exists());
    }

    #[test]
    fn hashing_matches_published_test_vectors() {
        let dir = tempfile::TempDir::new().unwrap();
        let a = dir.path().join("a");
        File::create(&a).unwrap().write_all(b"abc").unwrap();
        let (md5, sha) = hashes(&a).unwrap();
        assert_eq!(md5, "900150983cd24fb0d6963f7d28e17f72");
        assert_eq!(
            sha,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(md5_hex(b"abc"), md5);
    }
}
