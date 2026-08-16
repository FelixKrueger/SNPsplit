# SNPsplit fixture suite

A differential regression harness for `SNPsplit` and `tag2sort`. Every fixture runs the tool in a
private scratch directory and compares every file it produces against committed expected output.

Point `--impl` at any implementation to diff it against the same expected output, which makes
`fixtures/*/expected/` an executable specification rather than just a regression net.

## Running

```sh
test/run_tests.pl                 # every fixture
test/run_tests.pl se_basic hic    # named fixtures only
test/run_tests.pl --verbose       # echo each command
test/run_tests.pl --keep          # leave scratch dirs for inspection
test/run_tests.pl --impl /path/to/other/SNPsplit
```

Exit status is 0 when everything matches, 1 otherwise. `tag2sort` is taken from the same directory as
`--impl` and staged beside it, because `SNPsplit` invokes `$RealBin/tag2sort` with no override.

Requires `samtools` on `PATH` (or `--samtools PATH`) and core Perl only — no CPAN modules, matching
the tools under test.

## Re-recording the expected output

```sh
test/bin/make_fixtures.pl         # re-author reference, SNPs and SAM records
test/run_tests.pl --update        # re-record the expected output
git diff test/fixtures            # READ THIS
```

**Never run `--update` without reading the diff.** It will pin a bug just as happily as a fix. The
suite's value is entirely in the expected output being deliberate.

## Adding a fixture

Create `fixtures/<name>/` containing:

| File | Required | Purpose |
|---|---|---|
| `args` | yes | Extra command-line arguments, one line |
| `input.sam` | yes | Alignments; `input2.sam` adds a second input file |
| `snps.txt` | yes | Five-column SNP annotation |
| `chr1.fa` | yes | N-masked reference, used only at authoring time |
| `README` | yes | One line naming the code paths the fixture pins |
| `expected/` | generated | The expected output, recorded by --update |
| `feed_sam` | no | Marker: pass the input as `.sam` instead of converting to BAM |
| `poison_shim` | no | Marker: install a failing bare `samtools` first on `PATH` |
| `allow_empty` | no | Marker: this fixture legitimately produces no alignments |
| `exclude` | no | Glob per line: present in `manifest`, content not compared |

Prefer adding reads to `test/bin/make_fixtures.pl` rather than editing `input.sam` by hand, so MD tags
stay derived rather than asserted.

### Rules that are not optional

- **Use the registered option names.** They are `--SNP_file` and `--no_sorting`. `--snp_file` and
  `--no_sort` work only through Getopt::Long abbreviation, so relying on them makes every fixture
  depend on one Getopt::Long behaviour.
- **Every non-Hi-C fixture needs `--single_end` or `--paired`.** `check_for_bs` aborts the run
  otherwise, and it only reads `--single_end` inside its `@PG` branch, so `input.sam` must carry a
  `@PG` line. `--hic` bypasses the check entirely.
- **Never pass absolute paths.** The YAML records `infile` and `SNP_annotation` verbatim, so an
  absolute path bakes the random scratch directory into the recorded output.
- **Keep read counts small and avoid 32, 160 and 800.** `sprintf "%.2f"` half-way rounding can differ
  between glibc and libSystem when a percentage lands on exactly three decimals ending in 5.

## How MD tags are derived

Nothing in the toolchain validates an MD tag — `MD:Z:999N999` on a `10M` read round-trips through
samtools unchanged. Counter reconciliation cannot catch a tag that is self-consistent but points at
the wrong coordinate, because the counters reconcile perfectly around a wrong answer.

So `make_fixtures.pl` writes SAM records with correct POS/CIGAR/SEQ against a committed N-masked
reference and lets `samtools calmd` derive MD.

The reference and `snps.txt` can be checked against each other: every position masked to `N` in
`chr1.fa` is either listed in `snps.txt` or is one of the positions deliberately left out of it. Most
fixtures leave position 280 out — that is what gives `se_basic`'s `r_unlisted_n` read an `N` that
SNPsplit does not know about. `make_fixtures.pl` lists those positions under `unlisted`.

## Normalisation

The expected output is normalised text, not raw files:

- **BAM** → `samtools view -h`, with `@PG` `VN:`/`CL:` and `@HD` `SS:` values masked. `@PG` lines are
  kept because `check_for_bs` reads them to decide library type and bisulfite mode — they are a control
  input, not decoration. `SS:` is a recent sub-sort tag whose presence varies by samtools version.
- **YAML** → `version`, `date_run` and `command` masked. `command` is masked wholesale because a
  different implementation's argv differs legitimately.
- **Reports** → verbatim; they contain only basenames.
- **`run.log`** (stderr) → recorded and compared, with scratch and samtools paths masked. This is
  required, not extra: the hard-clipped-reads counter is reported *only* to stderr, so without it the
  suite cannot distinguish a hard-clip skip from an unmapped skip.
- **`manifest`** → sorted list of files produced. This is what pins `--conflicting`, `--singletons` and
  `--skip_tag2sort`, all of which gate whole files rather than content.
- **`exit_status`** → the run's exit code.

An unreadable or empty BAM normalises to `<unreadable>` rather than aborting the suite, so a fixture
can judge deliberately-broken output.

## Reading the reports

Two things look like bugs and are not:

- **Percentages can sum to more than 100%.** `$no_snp` is a subset of `$unassigned`, and
  `$non_N_containing` feeds `$unassigned` without feeding `$no_snp`.
- **Read counters do not sum to the total.** `$count` increments for every record including those later
  skipped. The identity that holds is:

  ```
  count - unmapped - hardclipped - noMD  ==  unassigned + genome1 + genome2 + conflicting
  ```

  `unmapped` is in the report, `hardclipped` is stderr-only, and reads with no `MD:Z` tag are counted
  nowhere at all.

## Fixtures that pin known defects

These assert current behaviour deliberately. Fixing any of them changes the expected output, which is a
real change to record, not a test failure.

| Fixture | Pinned behaviour | Issue |
|---|---|---|
| `hic_no_g1` | Hi-C with no G1/G1 pair prints a blank count and writes an empty YAML value | #90 |
| `output_dir` | `--output_dir` skips the sorting-report merge and leaves `SNPsplit_sort.yaml` behind | #91 |
| `multi_input` | `$snp_found`/`$no_snp_found` and `tag2sort`'s `%yaml` carry over between input files | #92 |
| `sam_output` | `--sam` produces no sorted output, because `tag2sort` only accepts `.bam` | #93 |
| `cigar_eqx` | `=`/`X` CIGAR operations abort the run; `H` is filtered but these are not | — |

`--sam` is broken end to end: `SNPsplit` writes `.allele_flagged.sam` and `tag2sort` accepts only
`.bam`. The `sam_output` fixture pins that, and is also the only fixture where the sorting stage fails,
so it is what tests SNPsplit aborting rather than reporting success. `sam_input` is a different
mechanism — `sam2bam_convert`, triggered by the input filename, not by `--sam`.

## Environment

`SNPSPLIT_NO_SLEEP=1` skips the tools' progress pauses. The harness sets it; it takes a single run from
roughly 2.2 s to 0.2 s and matters most for `--verbose`, which otherwise sleeps once per MD element.

The expected output in this repository was recorded with samtools 1.21 and perl 5.34.1 on macOS. CI runs
`ubuntu-24.04`. The same expected output has also been checked against samtools 1.22.1 with perl 5.32.1
and samtools 1.23.1 with perl 5.38.2 on Linux, and all three agree exactly — the samtools code that
writes SAM headers is unchanged from 1.19.1 to 1.21, and the volatile `@PG` and `@HD` fields are masked
anyway. If a future samtools does change what it writes, the fix is to re-record with `--update` and
check the diff is only the header line you expected.
