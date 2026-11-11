# The Rust suite

SNPsplit has been ported from Perl to Rust. The three tools keep their names, their options
and their output; what changes is that they are one binary with no Perl and no samtools
needed to run them.

!!! note "Status"

    The Rust build is released as SNPsplit 1.0, starting with 1.0.0-beta.1. The Perl version
    (0.9.x) is in maintenance freeze ([#167](https://github.com/FelixKrueger/SNPsplit/issues/167)).

## What is different

| | Perl | Rust |
|---|---|---|
| Runtime dependencies | Perl, samtools, gzip | none |
| Distribution | a folder of scripts | one binary, or a container image |
| Parallelism | none | `--parallel N`, opt-in |
| Fetching references | by hand | `--download`, opt-in |

Everything else is meant to be identical, and "meant to" is checked rather than asserted:
All committed fixtures but two produce byte-identical output, compared file by file after
the same normalisation the Perl suite uses. The two exceptions are described under
[Known differences](#known-differences).

## Installing

### Prebuilt binaries

Each Rust release, from 1.0.0-beta.1 on, attaches a `snpsplit-<version>-<platform>.tar.gz` to the
[Releases page](https://github.com/FelixKrueger/SNPsplit/releases) for Linux x86_64 and aarch64
(glibc 2.17 or newer) and macOS on Apple Silicon. Unpack it and put the folder on your `PATH`: it
holds `snpsplit` and the three classic names, which are links to it. On macOS, run
`xattr -d com.apple.quarantine snpsplit` once if the downloaded binary is blocked.

### Container image

```sh
docker run --rm --entrypoint SNPsplit ghcr.io/felixkrueger/snpsplit:beta --help
```

`:beta` is the newest prerelease, `:latest` the newest release, and each version has its own tag.
`SNPsplit`, `tag2sort` and `SNPsplit_genome_preparation` are all on the `PATH` inside, so a
pipeline that calls them by name is a drop-in. The image is distroless: no shell, no package
manager, no Perl.

### From source

```sh
git clone -b dev https://github.com/FelixKrueger/SNPsplit.git
cd SNPsplit
./rust/build-impl.sh
```

Requires a Rust toolchain (1.89 or newer) and a C compiler. This puts `SNPsplit`, `tag2sort` and
`SNPsplit_genome_preparation` in `rust/target/impl/`; add that folder to your `PATH`. The
`snpsplit` binary used in the examples below is `rust/target/release/snpsplit`. To build the image
yourself: `docker build -f rust/Dockerfile -t snpsplit .` from the same checkout.

## Invoking it

Either by the classic names, which is what existing pipelines do:

```sh
SNPsplit_genome_preparation --vcf_file mgp_REL2021_snps.vcf.gz --strain 129S1_SvImJ --reference_genome GRCm39/
SNPsplit --SNP_file all_SNPs_129S1_SvImJ_GRCm39.txt.gz sample.bam
```

or through the single binary:

```sh
snpsplit prepare --vcf_file ... --strain ...
snpsplit tag --SNP_file ... sample.bam
snpsplit sort sample.allele_flagged.bam
```

Both reach the same code. The classic names are links to one binary that decides which tool
to be from the name it was called by.

## `--parallel`

```sh
SNPsplit --SNP_file snps.txt --parallel 8 sample.bam
```

Default is 1, so an existing command line behaves exactly as it did. `--parallel 0` means
every available core, and `SNPSPLIT_PARALLEL` sets the default when the flag is absent, which
is convenient inside a container or under a scheduler.

**The output does not depend on the worker count.** Not "is deterministic for a given count",
which is weaker: the same input produces the same bytes at one worker and at sixteen. Both
fixture suites run at two different worker counts in CI against the same expected output.

How much it buys, measured on 2,000,000 reads:

| stage | 1 worker | 8 workers |
|---|---|---|
| Allele-tagging | 2.14s | 1.54s |
| Name sorting | 4.29s | 2.17s |

Worth being plain about: this is not a reason to switch. The per-read work is small next to
reading and writing the alignments, so most of what is left is I/O. The reason to switch is
not needing Perl and samtools installed.

## `--download`

```sh
SNPsplit_genome_preparation --download --strain 129S1_SvImJ
```

Fetches whatever of the two inputs is missing: the Mouse Genomes Project SNP VCF and a folder
of per-chromosome FastA files. `--download_dir PATH` says where (default
`./SNPsplit_references/`) and `--ensembl_release N` picks the genome release.

- Nothing contacts the network without the flag.
- A path you named with `--vcf_file` or `--reference_genome` always wins and is never
  downloaded over.
- Transfers resume rather than restart. The v8 VCF is 23 GB.
- The VCF is verified by decompressing it in full and requiring a VCF header, which catches a
  truncated transfer and an error page saved under a `.vcf.gz` name.
- URL, size and SHA-256 of everything fetched are written to `<download_dir>/manifest.txt`, so
  a later run can tell whether the remote file changed.

## Known differences

Two fixtures do not pass, and neither is a difference in what the tools compute.

`bam_write_fails` and `sam2bam_fails` make samtools fail to prove that a failed BAM write
aborts the run rather than leaving a truncated file behind. That property holds in the Rust
build, with the same messages, but there is no samtools subprocess for the fixtures to make
fail.

One thing the Rust build writes that the Perl does not: a `@PG` record naming itself in the
output header. The Perl records only the samtools invocations it shells out to, so the tool
that did the work never appears in its own output.

## Reporting a difference

If a run gives you something the Perl did not, that is a bug worth an issue, and the most
useful thing to include is the smallest input that shows it. The fixture suites are the
mechanism the port is held to, so a difference that can be turned into a fixture is a
difference that stays fixed.
