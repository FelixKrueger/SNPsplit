# SNPsplit Rust rewrite: design

Status: proposal. Nothing in this document is merged behaviour yet.

This is the design for porting SNPsplit from Perl to Rust, raised in
[#89](https://github.com/FelixKrueger/SNPsplit/issues/89). It exists so the shape of the
work is agreed before any of it lands, and so a reviewer can tell, at any point in the
stack, what is ported and what is not.

## Scope

All three tools:

| Tool | Perl lines | Fixture coverage today |
|---|---|---|
| `SNPsplit` | 1996 | `test/run_tests.pl`, 26 fixtures |
| `tag2sort` | 1311 | same suite (driven by `SNPsplit`) |
| `SNPsplit_genome_preparation` | 1788 | `test/run_genome_tests.pl`, 13 fixtures |

Issue #89 records a narrower option (`SNPsplit` + `tag2sort` as one binary, genome
preparation left in Perl) on the grounds that genome preparation is single-use and gains
nothing from being fast. That reasoning is about *performance*, and performance is not the
motivation here. The motivation is a single binary with no Perl and no samtools at runtime,
and a Perl genome preparation step would keep both dependencies alive for the one part of
the workflow that a user runs first. So: all three, one binary.

## What "done" means

Byte-identical output to Perl SNPsplit v0.9.0, proven by the existing fixture suites rather
than asserted.

Both runners already accept `--impl PATH`, and were written with this in mind ("Point
`--impl` at any implementation (Perl today, Rust later) to diff it against the same expected
output", `test/run_tests.pl:18`). They stage the implementation into a scratch directory,
run every fixture, normalise every file the run produces, and diff against committed
`expected/` trees. That is the acceptance gate for this port, unchanged.

Two details of the existing normalisation matter to the design:

- `test/run_tests.pl:384` rewrites `^(Samtools path:\s*)\S+` to `<samtools>`. The Rust
  implementation can print any single non-whitespace token there and still match.
- `test/run_genome_tests.pl:397` rewrites `<impl_name> line <digits>` to
  `<impl_name> line <n>`, masking the line number but keeping the rest of the Perl `die`
  suffix. A Rust error message has no such suffix at all, so every fixture whose `expected/`
  carries one has to be reconciled (see PR 3).

## Architecture

### Layout

```
rust/
  Cargo.toml          workspace
  Cargo.lock
  VERSION             suite version, single source of truth
  justfile            just impl, just test, just fixtures
  DESIGN.md           this file
  README.md           status journal: which module is ported, at what fidelity
  snpsplit/
    Cargo.toml
    src/
      main.rs         multicall dispatch
      lib.rs
      io/             noodles wrappers: SAM/BAM read, BAM write, name sort
      tag/            the SNPsplit tagging tool
      sort/           the tag2sort tool
      genome_prep/    the SNPsplit_genome_preparation tool
      report/         .txt report and .yaml report writers, shared
```

One crate with modules, not a crate per tool. The three tools share the SNP annotation
reader, the report writers and the BAM layer, and a workspace of three thin binary crates
around one fat library crate buys nothing but manifest churn. (Bismark's Rust suite started
as a crate per tool and converged on exactly this single-crate multicall shape.)

### Binary naming

One binary, `snpsplit`, dispatching on `argv[0]` and on a subcommand:

```
snpsplit tag ...     == SNPsplit ...
snpsplit sort ...    == tag2sort ...
snpsplit prepare ... == SNPsplit_genome_preparation ...
```

The three classic names are installed as symlinks to the one binary, so existing pipelines
and the fixture runners (which invoke `SNPsplit` and require `tag2sort` to sit beside it,
`test/run_tests.pl:42`) are drop-in.

During the port, `just impl` builds and links them into `rust/target/impl/`, so the gate is:

```
test/run_tests.pl --impl rust/target/impl/SNPsplit
test/run_genome_tests.pl --impl rust/target/impl/SNPsplit_genome_preparation
```

Whether the released binaries carry the classic names outright or a suffix during a beta
(as Bismark's `_rs` binaries did, so they could share a `PATH` with the Perl scripts) is a
packaging decision, deferred to PR 19.

### BAM I/O

Pure Rust via [noodles](https://github.com/zaeleus/noodles). No samtools subprocess, no
htslib build dependency.

This is the single largest behavioural change in the port, and it has consequences that are
not cosmetic:

1. **`--samtools_path` loses its job.** It stays accepted, so no existing command line
   breaks, but it is a no-op that prints a notice. The `Samtools path:` report line stays
   (masked by the runner anyway) so report layout is unchanged.
2. **Name sorting becomes ours.** `sort_by_name_paired_end` (`SNPsplit:1394`) shells out to
   `samtools sort -n`. In-memory sorting is not an option for real inputs, so PR 2 has to
   write an external merge sort with disk spill. This is the highest-risk item in the stack.
   Fallback, if it proves unreasonable: keep `samtools sort -n` as a subprocess for that one
   path, and record the exception in `rust/README.md` rather than hiding it.
3. **Four fixtures encode the samtools architecture** and cannot pass unmodified:
   `samtools_path`, `samtools_path_spaced` (both assert a poison shim is never reached),
   `bam_write_fails` and `sam2bam_fails` (both simulate a failing samtools wrapper). PR 18
   re-expresses each in terms the port can honour: the two path fixtures assert the no-op
   notice, and the two failure fixtures induce a real write failure (read-only output
   directory) instead of a fake samtools.

Every other output difference is a bug in the port, not a consequence of the design.

## The fidelity gate

A port that lands in twenty pieces spends most of its life partially complete, and the
honest question at every one of those points is "which fixtures pass?". Encoding the answer
in CI, rather than in a PR description, is what keeps the stack reviewable.

`test/rust_fixtures.txt` lists the fixtures the Rust implementation is expected to pass. CI
gains two jobs that run both suites against `rust/target/impl/` and fail on either side of
the line:

- a fixture on the list that fails, which is a regression;
- a fixture off the list that passes, which means the list is stale and someone forgot to
  claim their own work.

Failing in the second direction is the point. A plain skip list rots silently; this one
cannot. Each porting PR adds its fixtures to the list in the same commit that makes them
pass, so the diff always shows what was actually bought.

The Perl jobs stay exactly as they are and stay green throughout. The Perl implementation is
not touched by this stack.

## The stack

All PRs target the `rust` integration branch. Feature branches are `rs-*` rather than
`rust/*`, because git cannot hold both a `rust` ref and a `rust/...` ref at once.

The intent is that `rust` is merged to `master` once, as a single reviewed unit, with the
individual PRs serving as its readable history.

### Phase 0: foundation

| PR | Branch | Contents |
|---|---|---|
| 1 | `rs-scaffold` | Cargo workspace, `snpsplit` crate skeleton, multicall dispatch, `--version` and `--help` for all three names, `rust/README.md` journal, `just impl`, CI job for fmt/clippy/test |
| 2 | `rs-io` | noodles SAM/BAM reader and BAM writer, header and `@PG` handling, external name sort with disk spill, unit tests |
| 3 | `rs-harness` | `test/rust_fixtures.txt` plus the two CI jobs; reconcile the Perl `die`-shape masking in both runners. Perl suites stay green |

### Phase 1: genome preparation (13 fixtures, no BAM)

Ported first because it touches no alignments and has its own runner, so it proves the whole
chain (cargo, CI, allowlist, report fidelity) without waiting on `rs-io`.

| PR | Branch | Contents |
|---|---|---|
| 4 | `rs-gp-cli` | CLI surface, argument validation, `detect_chroms`, `detect_strains`; fixture `list_strains` |
| 5 | `rs-gp-vcf` | VCF filtering, high-confidence SNP selection, `SNPs_<strain>/chr*.txt`, gzipped `all_SNPs_*`, skipped-position reporting |
| 6 | `rs-gp-genome` | genome into memory, N-masking and full-sequence genomes; `nmask_basic`, `full_sequence`, `multi_chrom`, `multifasta`, `fasta_extension` |
| 7 | `rs-gp-dual` | dual hybrid path, `genotypes`, and the five failure fixtures (`chrom_name_mismatch`, `duplicate_chrom_name`, `empty_chromosome`, `empty_genome_folder`, `missing_genome`). Genome preparation suite complete at 13/13 |

### Phase 2: tag2sort

Ported before the tagger, because the Perl `SNPsplit` can drive a Rust `tag2sort`: the
runner stages both from `dirname(--impl)`, so a mixed directory is a valid gate.

| PR | Branch | Contents |
|---|---|---|
| 8 | `rs-sort-cli` | CLI, report and YAML writers |
| 9 | `rs-sort-se` | `process_single_end`; `se_basic`, `se_skipped` |
| 10 | `rs-sort-pe` | `princess_paired_end` and singleton handling; `pe_basic`, `pe_singletons`, `pe_namesort` |
| 11 | `rs-sort-hic` | `process_HiC_paired_end`; `hic`, `hic_no_g1` |

### Phase 3: SNPsplit tagging

| PR | Branch | Contents |
|---|---|---|
| 12 | `rs-tag-cli` | CLI, `read_snps` for both annotation formats including the per-chromosome format from #107, report and YAML |
| 13 | `rs-tag-se` | `process_masked_BAM_file` single-end, `score_snp_position`, `determine_N_position`, CIGAR handling; `cigar_complex`, `cigar_eqx`, `n_deletions` |
| 14 | `rs-tag-pe` | paired-end, the `--no_sorting` path; `multi_input`, `output_dir` |
| 15 | `rs-tag-hic` | Hi-C tagging |
| 16 | `rs-tag-bs` | `score_bisulfite_SNPs`, `check_for_bs` autodetection; `bisulfite_se`, `bismark_autodetect`, `bismark_autodetect_se` |
| 17 | `rs-tag-handoff` | in-process handoff from tagging to sorting, `--skip_tag2sort`, `sam_input`, `sam_output` |

### Phase 4: closing

| PR | Branch | Contents |
|---|---|---|
| 18 | `rs-fixtures-samtools` | re-express the four samtools-architecture fixtures. Alignment suite complete at 26/26 |
| 19 | `rs-packaging` | cargo metadata, release workflow, prebuilt binaries, container image |
| 20 | `rs-docs` | documentation site pages, README, CHANGELOG, and the body of the umbrella PR |

## Testing

Three layers, in increasing cost:

1. **Unit tests** inside each module, for the pieces with no I/O: CIGAR walking, SNP
   scoring, the bisulfite call, VCF field parsing. These are where a port actually goes
   wrong, and they are cheap enough to write per function.
2. **The two fixture suites** under `--impl`, gated by `test/rust_fixtures.txt`, run on every
   PR. This is the byte-identity contract.
3. **`test/bin/compare_versions.pl`**, the real-data harness, run by hand before the umbrella
   PR against a full mouse dataset. Fixtures are small by construction; a full run is the
   only thing that exercises memory behaviour and the external sort at scale.

## Error handling

Perl `die` messages are user-facing contract here: fixtures diff them. The port reproduces
the message text exactly, minus the ` at <script> line <n>.` suffix that Perl appends and
the runner masks. Exit statuses are diffed too (`expected/exit_status`), so they are
reproduced rather than approximated.

Internally, `anyhow` at the binary edge and typed errors inside modules. No `unwrap` on
anything derived from input: a malformed VCF line or a truncated BAM has to produce the
Perl message, not a panic.

## Deliberately out of scope

- Behaviour changes, new options, and output-format improvements. This is a port. Anything
  that would change output belongs in its own issue against the Perl version first, so the
  fixture suite records the change once and both implementations agree on it.
- Multithreading. The Perl tools are single-threaded, output order is part of the contract,
  and #89 is explicit that performance is not the motivation. Parallelism can come later,
  behind a flag, once byte-identity is established and can prove it did not break anything.
- Retiring the Perl scripts. That is Felix's call, after this lands, and it is a one-line
  change to make.
