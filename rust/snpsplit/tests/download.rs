//! The download path, against a local server.
//!
//! CI never contacts EBI or Ensembl. The server logs every request with its `Range` header, so
//! a test can tell a resume from a restart rather than only counting requests.

use std::collections::HashMap;
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::thread;

use snpsplit::genome_prep::download::{self, Manifest, Sources};
use tempfile::TempDir;

type Log = Arc<Mutex<Vec<(String, Option<u64>)>>>;

/// How the server answers for one path.
#[derive(Clone)]
enum Reply {
    File(Vec<u8>),
    /// Announces the whole file, sends half, then closes the connection.
    Drops(Vec<u8>),
    /// Ignores `Range` and always sends the whole file with a 200.
    IgnoresRange(Vec<u8>),
}

/// A Range-aware HTTP server for a fixed set of paths, on an ephemeral port.
struct Server {
    port: u16,
    log: Log,
}

impl Server {
    fn serving(files: Vec<(&str, Vec<u8>)>) -> Self {
        Self::replying(
            files
                .into_iter()
                .map(|(p, b)| (p.to_string(), Reply::File(b)))
                .collect(),
        )
    }

    fn replying(replies: Vec<(String, Reply)>) -> Self {
        let replies: Arc<HashMap<String, Reply>> = Arc::new(replies.into_iter().collect());
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind an ephemeral port");
        let port = listener.local_addr().unwrap().port();
        let log: Log = Arc::default();
        let thread_log = Arc::clone(&log);

        thread::spawn(move || {
            for stream in listener.incoming() {
                let Ok(stream) = stream else { break };
                handle(stream, &replies, &thread_log);
            }
        });

        Self { port, log }
    }

    fn url(&self, path: &str) -> String {
        format!("http://127.0.0.1:{}{path}", self.port)
    }

    fn requests(&self) -> Vec<(String, Option<u64>)> {
        self.log.lock().unwrap().clone()
    }
}

fn handle(mut stream: TcpStream, replies: &HashMap<String, Reply>, log: &Log) {
    let mut buffer = [0u8; 4096];
    let read = stream.read(&mut buffer).unwrap_or(0);
    let request = String::from_utf8_lossy(&buffer[..read]).to_string();
    let path = request.split_whitespace().nth(1).unwrap_or("/").to_string();
    let offset = request
        .lines()
        .find_map(|line| {
            let (name, value) = line.split_once(':')?;
            name.eq_ignore_ascii_case("range")
                .then(|| value.trim().strip_prefix("bytes="))
                .flatten()
        })
        .and_then(|value| value.trim_end_matches('-').parse::<u64>().ok());
    log.lock().unwrap().push((path.clone(), offset));

    let (body, drops, offset) = match replies.get(&path) {
        Some(Reply::File(b)) => (b, false, offset),
        Some(Reply::Drops(b)) => (b, true, offset),
        Some(Reply::IgnoresRange(b)) => (b, false, None),
        None => {
            let _ = stream.write_all(
                b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n",
            );
            return;
        }
    };

    // As Ensembl and EBI do: a range at or past the end is a 416 without Content-Range.
    let (status, slice) = match offset {
        Some(o) if o as usize >= body.len() => ("416 Range Not Satisfiable", &body[..0]),
        Some(o) => ("206 Partial Content", &body[o as usize..]),
        None => ("200 OK", &body[..]),
    };
    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
        slice.len()
    );
    let _ = stream.write_all(head.as_bytes());
    let sent = if drops {
        &slice[..slice.len() / 2]
    } else {
        slice
    };
    let _ = stream.write_all(sent);
    let _ = stream.flush();
}

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(data).unwrap();
    gz.finish().unwrap()
}

/// A VCF of two gzip members, so cutting it after the first leaves a file that decompresses
/// cleanly, as a BGZF file cut at a block boundary does.
fn two_member_vcf() -> (Vec<u8>, usize) {
    let mut first = String::from("##fileformat=VCFv4.2\n");
    for i in 0..500 {
        first.push_str(&format!("1\t{i}\t.\tA\tG\t.\tPASS\tAC=2\n"));
    }
    let mut second = String::new();
    for i in 500..1000 {
        second.push_str(&format!("1\t{i}\t.\tA\tG\t.\tPASS\tAC=2\n"));
    }
    let a = gzip(first.as_bytes());
    let cut = a.len();
    let mut body = a;
    body.extend(gzip(second.as_bytes()));
    (body, cut)
}

fn md5_line(data: &[u8]) -> Vec<u8> {
    format!("{}  mgp_REL2021_snps.vcf.gz\n", download::md5_hex(data)).into_bytes()
}

const VCF: &str = "/mgp_REL2021_snps.vcf.gz";
const MD5: &str = "/mgp_REL2021_snps.vcf.gz.md5";

fn vcf_server(body: &[u8]) -> Server {
    Server::serving(vec![(VCF, body.to_vec()), (MD5, md5_line(body))])
}

fn vcf_sources(server: &Server) -> Sources {
    Sources {
        vcf_url: server.url(VCF),
        vcf_md5_url: server.url(MD5),
        genome_base: server.url(""),
    }
}

fn part(path: &Path) -> PathBuf {
    PathBuf::from(format!("{}.part", path.display()))
}

fn ranges(server: &Server, path: &str) -> Vec<Option<u64>> {
    server
        .requests()
        .into_iter()
        .filter(|(p, _)| p == path)
        .map(|(_, r)| r)
        .collect()
}

/// `fetch_vcf` for a target the download chose, with a fresh manifest of `dir`.
fn fetch_vcf(dir: &Path, target: &Path, server: &Server) -> anyhow::Result<()> {
    let mut manifest = Manifest::load(dir).unwrap();
    download::fetch_vcf(
        target,
        false,
        &vcf_sources(server),
        &mut manifest,
        &mut Vec::new(),
    )
}

fn manifest_lines(dir: &Path) -> Vec<String> {
    std::fs::read_to_string(dir.join("manifest.txt"))
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect()
}

#[test]
fn a_fresh_fetch_lands_the_whole_file() {
    let server = Server::serving(vec![("/file.bin", b"0123456789".to_vec())]);
    let dir = TempDir::new().unwrap();
    let dest = dir.path().join("nested").join("file.bin");

    let bytes = download::fetch(&server.url("/file.bin"), &dest, &mut Vec::new()).unwrap();

    assert_eq!(bytes, 10);
    assert_eq!(std::fs::read(&dest).unwrap(), b"0123456789");
}

#[test]
fn a_missing_file_fails_loudly_and_names_the_url() {
    let server = Server::serving(vec![]);
    let dir = TempDir::new().unwrap();
    let url = server.url("/missing");
    let Err(e) = download::fetch(&url, &dir.path().join("f"), &mut Vec::new()) else {
        panic!("a 404 should not be treated as a download");
    };
    assert!(
        e.to_string().contains(&url),
        "error did not name the URL: {e}"
    );
}

#[test]
fn a_range_past_the_end_counts_as_transferred() {
    let server = Server::serving(vec![("/file.bin", b"0123456789".to_vec())]);
    let dir = TempDir::new().unwrap();
    let dest = dir.path().join("file.bin");
    std::fs::write(&dest, b"0123456789").unwrap();

    let bytes = download::fetch(&server.url("/file.bin"), &dest, &mut Vec::new()).unwrap();

    assert_eq!(bytes, 10);
    assert_eq!(ranges(&server, "/file.bin"), vec![Some(10)]);
}

#[test]
fn a_dropped_transfer_keeps_the_part_for_the_next_run() {
    let (body, _) = two_member_vcf();
    let server = Server::replying(vec![
        (VCF.to_string(), Reply::Drops(body.clone())),
        (MD5.to_string(), Reply::File(md5_line(&body))),
    ]);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");

    let e = fetch_vcf(dir.path(), &target, &server).unwrap_err();

    assert!(
        format!("{e:#}").contains("run the same command again"),
        "{e:#}"
    );
    assert!(!target.exists());
    assert_eq!(
        std::fs::read(part(&target)).unwrap(),
        &body[..body.len() / 2]
    );
    assert!(manifest_lines(dir.path()).is_empty());

    let rerun = vcf_server(&body);
    fetch_vcf(dir.path(), &target, &rerun).unwrap();
    assert_eq!(ranges(&rerun, VCF), vec![Some((body.len() / 2) as u64)]);
    assert_eq!(std::fs::read(&target).unwrap(), body);
}

#[test]
fn an_interrupted_vcf_is_resumed_and_checked() {
    let (body, _) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    let half = body.len() / 2;
    std::fs::write(part(&target), &body[..half]).unwrap();

    fetch_vcf(dir.path(), &target, &server).unwrap();

    assert_eq!(std::fs::read(&target).unwrap(), body);
    assert!(!part(&target).exists());
    assert_eq!(ranges(&server, VCF), vec![Some(half as u64)]);
    let lines = manifest_lines(dir.path());
    assert_eq!(lines.len(), 1);
    assert!(
        lines[0].starts_with('/'),
        "paths are recorded absolute: {}",
        lines[0]
    );
    assert!(Manifest::load(dir.path()).unwrap().lists(&target));
}

#[test]
fn a_complete_part_gets_a_416_and_is_then_checked() {
    let (body, _) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    std::fs::write(part(&target), &body).unwrap();

    fetch_vcf(dir.path(), &target, &server).unwrap();

    assert_eq!(ranges(&server, VCF), vec![Some(body.len() as u64)]);
    assert_eq!(std::fs::read(&target).unwrap(), body);
}

#[test]
fn a_vcf_left_truncated_by_1_0_0_is_resumed() {
    let (body, cut) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    // Cut at a member boundary: it decompresses cleanly and has a header, but is incomplete.
    std::fs::write(&target, &body[..cut]).unwrap();

    fetch_vcf(dir.path(), &target, &server).unwrap();

    assert_eq!(std::fs::read(&target).unwrap(), body);
    assert_eq!(ranges(&server, VCF), vec![Some(cut as u64)]);
}

#[test]
fn a_complete_vcf_from_1_0_0_is_checked_not_refetched() {
    let (body, _) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    std::fs::write(&target, &body).unwrap();

    fetch_vcf(dir.path(), &target, &server).unwrap();

    assert!(ranges(&server, VCF).is_empty(), "the VCF was fetched again");
    assert_eq!(ranges(&server, MD5).len(), 1);
    assert!(Manifest::load(dir.path()).unwrap().lists(&target));
}

#[test]
fn a_listed_vcf_needs_no_request_unless_its_size_changed() {
    let (body, _) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    fetch_vcf(dir.path(), &target, &server).unwrap();
    let before = server.requests().len();

    fetch_vcf(dir.path(), &target, &server).unwrap();
    assert_eq!(server.requests().len(), before);

    let mut grown = body.clone();
    grown.extend(b"garbage");
    std::fs::write(&target, &grown).unwrap();
    fetch_vcf(dir.path(), &target, &server).unwrap();
    assert_eq!(std::fs::read(&target).unwrap(), body);
}

#[test]
fn an_error_page_left_as_part_is_restarted_not_resumed() {
    let (body, _) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    std::fs::write(part(&target), b"<html>503 Service Unavailable</html>").unwrap();

    fetch_vcf(dir.path(), &target, &server).unwrap();

    assert_eq!(std::fs::read(&target).unwrap(), body);
    assert_eq!(ranges(&server, VCF), vec![None]);
}

#[test]
fn a_resume_that_does_not_check_out_is_restarted_once() {
    let (body, _) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    // A gzip file that is not the start of this VCF: resuming it cannot give the right MD5.
    std::fs::write(
        part(&target),
        gzip(b"##fileformat=VCFv4.2\nsomething else\n"),
    )
    .unwrap();

    fetch_vcf(dir.path(), &target, &server).unwrap();

    assert_eq!(std::fs::read(&target).unwrap(), body);
    let r = ranges(&server, VCF);
    assert_eq!(r.len(), 2);
    assert!(
        r[0].is_some() && r[1].is_none(),
        "expected a resume, then a restart: {r:?}"
    );
}

#[test]
fn a_vcf_that_never_checks_out_fails_and_leaves_nothing() {
    let (body, _) = two_member_vcf();
    let server = Server::serving(vec![
        (VCF, b"<html>not a VCF</html>".to_vec()),
        (MD5, md5_line(&body)),
    ]);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");

    let e = fetch_vcf(dir.path(), &target, &server).unwrap_err();

    assert!(e.to_string().contains(&server.url(VCF)), "{e}");
    assert!(!target.exists() && !part(&target).exists());
}

#[test]
fn a_missing_or_garbled_md5_stops_before_the_vcf_is_fetched() {
    let (body, _) = two_member_vcf();
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");

    let no_md5 = Server::serving(vec![(VCF, body.clone())]);
    let e = fetch_vcf(dir.path(), &target, &no_md5).unwrap_err();
    assert!(e.to_string().contains(&no_md5.url(MD5)), "{e}");

    let html_md5 = Server::serving(vec![(VCF, body.clone()), (MD5, b"<html></html>".to_vec())]);
    let e = fetch_vcf(dir.path(), &target, &html_md5).unwrap_err();
    assert!(e.to_string().contains("does not hold an MD5"), "{e}");

    assert!(ranges(&no_md5, VCF).is_empty() && ranges(&html_md5, VCF).is_empty());
}

#[test]
fn a_named_vcf_that_fails_the_check_is_left_unchanged() {
    let (body, _) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let mine = dir.path().join("mine");
    std::fs::create_dir(&mine).unwrap();
    let named = mine.join("mgp_REL2021_snps.vcf.gz");
    let own = gzip(b"##fileformat=VCFv4.2\nmy own edit\n");
    std::fs::write(&named, &own).unwrap();
    let mut manifest = Manifest::load(dir.path()).unwrap();

    let e = download::fetch_vcf(
        &named,
        true,
        &vcf_sources(&server),
        &mut manifest,
        &mut Vec::new(),
    )
    .unwrap_err();

    assert!(e.to_string().contains("left unchanged"), "{e}");
    assert_eq!(std::fs::read(&named).unwrap(), own);
    assert!(!part(&named).exists());
    assert!(ranges(&server, VCF).is_empty());
}

#[test]
fn a_symlinked_vcf_that_fails_the_check_is_not_written_through() {
    let (body, cut) = two_member_vcf();
    let server = vcf_server(&body);
    let dir = TempDir::new().unwrap();
    let real = dir.path().join("elsewhere.vcf.gz");
    std::fs::write(&real, &body[..cut]).unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    std::os::unix::fs::symlink(&real, &target).unwrap();

    assert!(fetch_vcf(dir.path(), &target, &server).is_err());
    assert_eq!(std::fs::read(&real).unwrap(), &body[..cut]);
}

fn chromosome_files(build: &str, names: &[&str]) -> Vec<(String, Vec<u8>)> {
    names
        .iter()
        .map(|c| {
            (
                format!("/Mus_musculus.{build}.dna.chromosome.{c}.fa.gz"),
                gzip(format!(">{c}\n{}\n", "ACGT".repeat(500)).as_bytes()),
            )
        })
        .collect()
}

fn genome_server(names: &[&str]) -> Server {
    let files = chromosome_files("GRCm39", names);
    Server::serving(files.iter().map(|(p, b)| (p.as_str(), b.clone())).collect())
}

fn chroms(names: &[&str]) -> Vec<String> {
    names.iter().map(|s| s.to_string()).collect()
}

fn requested_chromosomes(server: &Server) -> Vec<String> {
    server
        .requests()
        .into_iter()
        .map(|(p, _)| {
            p.trim_start_matches("/Mus_musculus.GRCm39.dna.chromosome.")
                .trim_end_matches(".fa.gz")
                .to_string()
        })
        .collect()
}

fn fetch_genome(dir: &Path, genome: &Path, server: &Server, names: &[&str]) -> anyhow::Result<()> {
    let mut manifest = Manifest::load(dir).unwrap();
    download::fetch_genome(
        genome,
        &server.url(""),
        "GRCm39",
        115,
        &chroms(names),
        &mut manifest,
        &mut Vec::new(),
    )
}

fn unpacked(c: &str) -> String {
    format!(">{c}\n{}\n", "ACGT".repeat(500))
}

#[test]
fn an_interrupted_genome_fetches_only_what_is_missing() {
    let all = ["1", "2", "3", "4", "5", "6"];
    let server = genome_server(&all);
    let dir = TempDir::new().unwrap();
    let genome = dir.path().join("GRCm39");
    std::fs::create_dir_all(&genome).unwrap();
    for c in &all[..4] {
        std::fs::write(genome.join(format!("{c}.fa")), unpacked(c)).unwrap();
    }

    fetch_genome(dir.path(), &genome, &server, &all).unwrap();

    assert_eq!(requested_chromosomes(&server), vec!["5", "6"]);
    for c in &all {
        assert_eq!(
            std::fs::read_to_string(genome.join(format!("{c}.fa"))).unwrap(),
            unpacked(c)
        );
    }
    assert_eq!(manifest_lines(dir.path()).len(), 2);
    assert_eq!(
        std::fs::read_to_string(genome.join(download::MARKER)).unwrap(),
        "GRCm39 115\n"
    );
}

#[test]
fn a_partial_archive_is_resumed() {
    let all = ["1"];
    let server = genome_server(&all);
    let files = chromosome_files("GRCm39", &all);
    let dir = TempDir::new().unwrap();
    let genome = dir.path().join("GRCm39");
    std::fs::create_dir_all(&genome).unwrap();
    let half = files[0].1.len() / 2;
    let gz_part = genome.join(format!("{}.part", files[0].0.trim_start_matches('/')));
    std::fs::write(&gz_part, &files[0].1[..half]).unwrap();

    fetch_genome(dir.path(), &genome, &server, &all).unwrap();

    assert_eq!(ranges(&server, &files[0].0), vec![Some(half as u64)]);
    assert_eq!(
        std::fs::read_to_string(genome.join("1.fa")).unwrap(),
        unpacked("1")
    );
    assert!(!gz_part.exists());
}

#[test]
fn a_partial_archive_that_is_not_gzip_is_restarted() {
    let all = ["1"];
    let server = genome_server(&all);
    let files = chromosome_files("GRCm39", &all);
    let dir = TempDir::new().unwrap();
    let genome = dir.path().join("GRCm39");
    std::fs::create_dir_all(&genome).unwrap();
    let gz_part = genome.join(format!("{}.part", files[0].0.trim_start_matches('/')));
    std::fs::write(&gz_part, b"<html>busy</html>").unwrap();

    fetch_genome(dir.path(), &genome, &server, &all).unwrap();

    assert_eq!(ranges(&server, &files[0].0), vec![None]);
    assert_eq!(
        std::fs::read_to_string(genome.join("1.fa")).unwrap(),
        unpacked("1")
    );
}

#[test]
fn a_resumed_archive_that_does_not_unpack_is_fetched_again() {
    let all = ["1"];
    let server = genome_server(&all);
    let files = chromosome_files("GRCm39", &all);
    let dir = TempDir::new().unwrap();
    let genome = dir.path().join("GRCm39");
    std::fs::create_dir_all(&genome).unwrap();
    // A gzip start from some other file: the resumed archive cannot unpack.
    let other = gzip(b">X\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n");
    let gz_part = genome.join(format!("{}.part", files[0].0.trim_start_matches('/')));
    std::fs::write(&gz_part, &other[..other.len() / 2]).unwrap();

    fetch_genome(dir.path(), &genome, &server, &all).unwrap();

    let r = ranges(&server, &files[0].0);
    assert_eq!(r.len(), 2, "{r:?}");
    assert!(
        r[0].is_some() && r[1].is_none(),
        "expected a resume, then a restart: {r:?}"
    );
    assert_eq!(
        std::fs::read_to_string(genome.join("1.fa")).unwrap(),
        unpacked("1")
    );
}

#[test]
fn a_chromosome_whose_archive_is_still_there_is_unpacked_again() {
    let all = ["1", "2"];
    let server = genome_server(&all);
    let dir = TempDir::new().unwrap();
    let genome = dir.path().join("GRCm39");
    std::fs::create_dir_all(&genome).unwrap();
    let files = chromosome_files("GRCm39", &all);
    // As 1.0.0 left it when stopped while unpacking chromosome 2.
    std::fs::write(genome.join("1.fa"), unpacked("1")).unwrap();
    std::fs::write(genome.join("2.fa"), ">2\nAC").unwrap();
    std::fs::write(genome.join(files[1].0.trim_start_matches('/')), &files[1].1).unwrap();

    fetch_genome(dir.path(), &genome, &server, &all).unwrap();

    assert_eq!(
        std::fs::read_to_string(genome.join("2.fa")).unwrap(),
        unpacked("2")
    );
    assert_eq!(requested_chromosomes(&server), vec!["2"]);
    assert_eq!(
        ranges(&server, &files[1].0),
        vec![Some(files[1].1.len() as u64)]
    );
}

#[test]
fn a_missing_chromosome_stops_the_run_and_a_rerun_finishes_it() {
    let all = ["1", "2", "3"];
    let dir = TempDir::new().unwrap();
    let genome = dir.path().join("GRCm39");

    let first = genome_server(&all[..2]);
    let e = fetch_genome(dir.path(), &genome, &first, &all).unwrap_err();
    assert!(e.to_string().contains("chromosome.3.fa.gz"), "{e}");
    let listed = Manifest::load(dir.path()).unwrap();
    assert!(listed.contains(&genome.join("1.fa")) && listed.contains(&genome.join("2.fa")));
    assert_eq!(
        download::genome_owner(&genome, false, "GRCm39"),
        download::Owner::Ours
    );

    let second = genome_server(&all);
    fetch_genome(dir.path(), &genome, &second, &all).unwrap();
    assert_eq!(requested_chromosomes(&second), vec!["3"]);
    assert_eq!(manifest_lines(dir.path()).len(), 3);
}

#[test]
fn a_resume_answered_with_the_whole_file_starts_again() {
    let (body, _) = two_member_vcf();
    let server = Server::replying(vec![
        (VCF.to_string(), Reply::IgnoresRange(body.clone())),
        (MD5.to_string(), Reply::File(md5_line(&body))),
    ]);
    let dir = TempDir::new().unwrap();
    let target = dir.path().join("mgp_REL2021_snps.vcf.gz");
    std::fs::write(part(&target), &body[..body.len() / 2]).unwrap();

    fetch_vcf(dir.path(), &target, &server).unwrap();

    assert_eq!(
        std::fs::read(&target).unwrap(),
        body,
        "the part was appended to"
    );
}

#[test]
fn a_complete_folder_is_used_without_writing_to_it() {
    let all = ["1", "2"];
    let server = genome_server(&all);
    let dir = TempDir::new().unwrap();
    let genome = dir.path().join("GRCm39");
    std::fs::create_dir_all(&genome).unwrap();
    for c in &all {
        std::fs::write(genome.join(format!("{c}.fa")), unpacked(c)).unwrap();
    }
    std::fs::write(genome.join(download::MARKER), "GRCm39 114\n").unwrap();
    let before = std::fs::read_dir(&genome).unwrap().count();
    let mut readonly = std::fs::metadata(&genome).unwrap().permissions();
    readonly.set_readonly(true);
    std::fs::set_permissions(&genome, readonly.clone()).unwrap();

    // Release 115 against a complete folder from 114: nothing is missing, so nothing is mixed.
    let result = fetch_genome(dir.path(), &genome, &server, &all);

    #[allow(clippy::permissions_set_readonly_false)]
    readonly.set_readonly(false);
    std::fs::set_permissions(&genome, readonly).unwrap();
    result.unwrap();
    assert!(server.requests().is_empty());
    assert_eq!(std::fs::read_dir(&genome).unwrap().count(), before);
    assert_eq!(
        std::fs::read_to_string(genome.join(download::MARKER)).unwrap(),
        "GRCm39 114\n"
    );
}
