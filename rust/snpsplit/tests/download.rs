//! The download path, against a local server.
//!
//! CI never contacts EBI or Ensembl. What is worth testing here is the transfer itself:
//! whether a fresh fetch lands, whether a partial file is resumed rather than restarted, and
//! whether a missing file fails loudly. The URLs the tool uses in anger are covered by unit
//! tests of the source table.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::Path;
use std::thread;

use snpsplit::genome_prep::download;
use tempfile::TempDir;

/// A one-file HTTP server that understands `Range`, on an ephemeral port.
struct Server {
    port: u16,
}

impl Server {
    fn serving(body: &'static [u8]) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind an ephemeral port");
        let port = listener.local_addr().unwrap().port();

        thread::spawn(move || {
            for stream in listener.incoming() {
                let Ok(stream) = stream else { break };
                handle(stream, body);
            }
        });

        Self { port }
    }

    fn url(&self, path: &str) -> String {
        format!("http://127.0.0.1:{}{path}", self.port)
    }
}

fn handle(mut stream: TcpStream, body: &[u8]) {
    let mut buffer = [0u8; 2048];
    let read = stream.read(&mut buffer).unwrap_or(0);
    let request = String::from_utf8_lossy(&buffer[..read]).to_string();

    if request.contains("/missing") {
        let _ = stream
            .write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n");
        return;
    }

    // A Range header gets the tail and a 206, which is what resuming depends on.
    let offset = request
        .lines()
        .find_map(|line| line.strip_prefix("Range: bytes="))
        .and_then(|value| value.trim_end_matches('-').parse::<usize>().ok());

    let (status, slice) = match offset {
        Some(offset) if offset < body.len() => ("206 Partial Content", &body[offset..]),
        Some(_) => ("206 Partial Content", &body[body.len()..]),
        None => ("200 OK", body),
    };

    let head = format!(
        "HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
        slice.len()
    );
    let _ = stream.write_all(head.as_bytes());
    let _ = stream.write_all(slice);
    let _ = stream.flush();
}

const BODY: &[u8] = b"0123456789abcdefghijklmnopqrstuvwxyz";

#[test]
fn a_fresh_fetch_lands_the_whole_file() {
    let server = Server::serving(BODY);
    let dir = TempDir::new().unwrap();
    let dest = dir.path().join("nested").join("file.bin");

    let bytes = download::fetch(&server.url("/file.bin"), &dest, &mut Vec::new()).unwrap();

    assert_eq!(bytes, BODY.len() as u64);
    assert_eq!(std::fs::read(&dest).unwrap(), BODY);
}

/// The reason resuming exists: the MGP VCF is large enough that losing a transfer to a
/// dropped connection is a real cost, and a half-file from a previous run is the case worth
/// handling.
#[test]
fn a_partial_file_is_resumed_rather_than_restarted() {
    let server = Server::serving(BODY);
    let dir = TempDir::new().unwrap();
    let dest = dir.path().join("file.bin");

    // Ten bytes already on disk, as a dropped transfer would leave.
    std::fs::write(&dest, &BODY[..10]).unwrap();

    let mut log = Vec::new();
    let bytes = download::fetch(&server.url("/file.bin"), &dest, &mut log).unwrap();

    assert_eq!(bytes, BODY.len() as u64);
    assert_eq!(std::fs::read(&dest).unwrap(), BODY);
    assert!(
        String::from_utf8_lossy(&log).contains("Resuming"),
        "a resumed transfer should say so: {}",
        String::from_utf8_lossy(&log)
    );
}

#[test]
fn a_missing_file_fails_loudly_and_names_the_url() {
    let server = Server::serving(BODY);
    let dir = TempDir::new().unwrap();
    let dest = dir.path().join("file.bin");

    let url = server.url("/missing");
    let Err(e) = download::fetch(&url, &dest, &mut Vec::new()) else {
        panic!("a 404 should not be treated as a download");
    };
    assert!(
        e.to_string().contains(&url),
        "error did not name the URL: {e}"
    );
}

/// A truncated transfer has to be caught before the file is used, not during the run that
/// uses it.
#[test]
fn a_truncated_vcf_does_not_verify() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("truncated.vcf.gz");

    let file = std::fs::File::create(&path).unwrap();
    let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    gz.write_all(b"##fileformat=VCFv4.2\n").unwrap();
    for i in 0..1000 {
        writeln!(gz, "1\t{i}\t.\tA\tG\t.\tPASS\tAC=2").unwrap();
    }
    gz.finish().unwrap();

    // Cut the tail off, as an interrupted transfer would.
    let full = std::fs::read(&path).unwrap();
    std::fs::write(&path, &full[..full.len() / 2]).unwrap();

    assert!(
        download::verify_vcf(Path::new(&path)).is_err(),
        "a truncated archive should not verify"
    );
}
