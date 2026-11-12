[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/snpsplit/README.html)

<p align="center"> <img title="SNPsplit" id="logo_img" src="https://raw.githubusercontent.com/FelixKrueger/SNPsplit/master/docs/images/SNPsplit.png" width=300></p>

## Allele-specific alignment sorting

SNPsplit sorts reads from hybrid samples, such as mouse crosses, by parental origin. Reads are
aligned to a genome in which known SNP positions are masked with N, and SNPsplit then assigns each
read that covers a SNP to genome 1 or genome 2. It works for most library types, including
RNA-seq, ChIP-seq, ATAC-seq, bisulfite-seq (Bismark) and Hi-C (HiCUP).

> **Documentation: https://felixkrueger.github.io/SNPsplit**

SNPsplit 1.0 is a rewrite in Rust: one binary that gives the same output as the Perl version, is
3 to 5 times faster and uses about a tenth of the memory, and needs nothing else installed.

## Installation

**Prebuilt binaries.** Download the archive for your platform from the
[Releases page](https://github.com/FelixKrueger/SNPsplit/releases), unpack it and put the folder on
your `PATH`:

| Platform | Archive |
|---|---|
| Linux x86_64 (glibc 2.17 or newer) | `snpsplit-<version>-linux-x86_64.tar.gz` |
| Linux aarch64 (glibc 2.17 or newer) | `snpsplit-<version>-linux-aarch64.tar.gz` |
| macOS Apple Silicon | `snpsplit-<version>-macos-aarch64.tar.gz` |

It holds `SNPsplit`, `tag2sort` and `SNPsplit_genome_preparation`. On macOS, run
`xattr -d com.apple.quarantine SNPsplit` once if the downloaded binary is blocked.

**Container image:**

```sh
docker run --rm --entrypoint SNPsplit ghcr.io/felixkrueger/snpsplit:latest --help
```

**crates.io:** `cargo install snpsplit` installs a single `snpsplit` binary, run by subcommand
(`snpsplit prepare`, `snpsplit tag`, `snpsplit sort`). It needs a Rust toolchain (1.89 or newer)
and a C compiler.

**Bioconda:** `conda install -c bioconda snpsplit` currently installs the Perl version; the recipe
moves to 1.0 shortly after its release.

**From source:** see [the Rust docs](https://felixkrueger.github.io/SNPsplit/rust/).

## Quick start

```sh
# 1. Build an N-masked genome and the SNP file for one strain (Mouse Genomes Project VCF)
SNPsplit_genome_preparation --vcf_file mgp_REL2021_snps.vcf.gz --strain CAST_EiJ --reference_genome GRCm39/

# 2. Index CAST_EiJ_N-masked/ with your aligner and align the reads to it, then sort the alignments
SNPsplit --snp_file all_SNPs_CAST_EiJ_GRCm39.txt.gz sample.bam
```

This writes `sample.genome1.bam`, `sample.genome2.bam`, `sample.unassigned.bam` and a report.
Use `--paired`, `--bisulfite` or `--hic` to match the data, and `--parallel N` to use N workers.
`SNPsplit_genome_preparation --download` fetches the VCF and the reference if they are missing.
The documentation covers dual hybrid genomes, every option and the reports.

## The Perl version

The Perl version (0.9.x) is in maintenance freeze and gets critical fixes only
([#167](https://github.com/FelixKrueger/SNPsplit/issues/167)). To use it, download the
[0.9.2 release](https://github.com/FelixKrueger/SNPsplit/releases/tag/0.9.2); it needs
[Samtools](http://samtools.sourceforge.net/) on the `PATH`. In this repository the Perl scripts are
in `legacy_perl/` and the Rust is in `rust/`.

## Citation

Krueger F and Andrews SR. SNPsplit: Allele-specific splitting of alignments between genomes with
known SNP genotypes. *F1000Research* 2016, 5:1479. https://f1000research.com/articles/5-1479/v2

## Credits

SNPsplit was written by Felix Krueger at [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk),
now part of Altos Bioinformatics. The Rust port was started by Benjamin Demaille
([#89](https://github.com/FelixKrueger/SNPsplit/issues/89)). The original project page is at the
[Babraham Institute](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/).

<p align="center"> <img title="Babraham Bioinformatics" id="logo_img" src="https://raw.githubusercontent.com/FelixKrueger/SNPsplit/master/docs/images/bioinformatics_logo.png" width=200></p>
