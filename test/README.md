# SNPsplit fixture suite

Differential regression harnesses for all three tools. Every fixture runs the tool in a private
scratch directory and compares every file it produces against committed expected output.

Point `--impl` at any implementation to diff it against the same expected output, which makes
`*/expected/` an executable specification rather than just a regression net.

| Runner | Tools | Fixtures | Needs |
|---|---|---|---|
| `run_tests.pl` | `SNPsplit`, `tag2sort` | `fixtures/`, 26 | `samtools`, `gzip` |
| `run_genome_tests.pl` | `SNPsplit_genome_preparation` | `genome_fixtures/`, 31 | `gzip` only |

Two runners rather than one: the tools disagree on nearly every tool-specific decision a runner makes.
Genome preparation needs no samtools, takes a VCF plus a reference *directory* instead of alignments,
has no companion binary to stage, and needs three normalisations the other suite has no use for. The
tool-agnostic halves of the two runners are byte-identical so they can be de-duplicated into a shared
module mechanically; the four subs that could not be copied verbatim are named at their definitions
with the reason.

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
| `break_tag2sort` | no | Marker: stage a sorter that always fails, to test the abort |
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
reference and lets `samtools calmd` derive MD. The runner re-derives them again on every run and dies on
any mismatch, so a hand-edited `input.sam` is caught rather than quietly recorded.

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
- **`exit_status`** → `0` or `nonzero`, not the number. Perl's `die` exits with `$!`, which reflects the
  last failed syscall and so changes depending on what ran before; the die message itself is compared as
  part of `run.log`, so the number carried no information and made fixtures intermittent.

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

## Fixtures worth knowing about

**`cigar_eqx`** asserts that `=` and `X` CIGAR operations abort the run. `H` is filtered out earlier but
these are not, so output from an aligner run with `--eqx` cannot be processed. Fixing that would change
the expected output, which is a real change to record rather than a test failure.

**`tag2sort_fails`** replaces the staged sorter with one that always exits non-zero, so it tests that
SNPsplit aborts when the sorting stage fails, whatever the reason. It is written that way deliberately:
that behaviour used to be covered only as a side effect of `--sam` being broken, and the coverage
disappeared the moment `--sam` was fixed. Coverage that depends on a bug still being present is coverage
on loan.

**`sam_input` and `sam_output`** cover two different mechanisms. `sam_input` exercises
`sam2bam_convert`, triggered by the input *filename*; `sam_output` exercises `--sam`, which sets the
output format.

`hic_no_g1`, `output_dir`, `multi_input` and `sam_output` used to pin defects and now assert correct
behaviour, since #90 to #93 are fixed.

**`gp_generated`** takes its `chr1.fa` and `snps.txt` from a real `SNPsplit_genome_preparation` run at
authoring time instead of writing them by hand, so the file one tool writes and the other reads is
checked rather than assumed. Same reasoning as the MD tags. Regenerate it with
`test/bin/make_fixtures.pl`, which runs the genome preparation script itself.

The annotation comes from `all_SNPs_<strain>_<build>.txt.gz`, not from `SNPs_<strain>/chr*.txt`: the
per-chromosome file opens with a `>chr` header line and `read_snps` has no header skip, so reading it
yields uninitialised-value warnings and a SNP stored under an empty chromosome name. It is sorted by
position, because the archive is written from a hash and would otherwise land in a different order every
time anyone regenerates the fixture.

It carries its own reference rather than sharing `%STD`: genome preparation masks exactly what it
filters, so it cannot produce the unlisted `N` at position 280 that `se_basic` depends on.

# Genome preparation suite

`run_genome_tests.pl` drives `SNPsplit_genome_preparation` over `genome_fixtures/`. Same flags as the
other runner minus `--samtools`, and it injects no arguments: each fixture's `args` file is the whole
command line, because this tool takes three independent mandatory-ish inputs rather than one.

```sh
test/bin/make_genome_fixtures.pl   # re-author references, VCFs and SNP files
test/run_genome_tests.pl --update  # re-record
git diff test/genome_fixtures      # READ THIS
```

`make_genome_fixtures.pl` rewrites only the inputs it owns and never touches `expected/`.

## Fixture layout

| File | Required | Purpose |
|---|---|---|
| `args` | yes | The complete argument list, one line |
| `README` | yes | One line naming the code paths pinned |
| `genome/` | no | Reference FastA files, copied verbatim. Omitted by fixtures that abort before reading it |
| `<name>.vcf` | no | Staged under its own basename, so `vcf_v7` can be called `mgp_REL2005_snps_indels.vcf` |
| `snps_in/` | no | Copied to the scratch root: a pre-made `SNPs_<strain>/` tree for `--skip_filtering`. A file here named `*.gz` holds **plain text** and is compressed on the way in, so no binary lands in the repository |
| `pre_existing` | no | Filename per line, created empty before the run |
| `gzip_vcf` | no | Marker: stage the VCF gzipped |
| `poison_gzip` | no | Marker: install a failing `gzip` first on `PATH` |
| `allow_empty` | no | Marker: this fixture legitimately generates no genome |
| `exclude` | no | Glob per line: in `manifest`, content not compared. **Cannot reach `stdout.log`** |

## Rules that are not optional

- **Never type a VCF `REF` base.** Read it out of the reference the VCF is written against.
  `create_modified_chromosome` skips a SNP whose `REF` disagrees, and the counter that records the skip
  is never reported, so a mistyped base produces a fixture that records a plausible zero-SNP result and
  looks correct. `make_genome_fixtures.pl` asserts every `REF` after writing it.
- **Give every hash-ordered output at least five records.** Below that the order cannot vary, or repeats
  by chance, and the fixture passes without exercising the normalisation it exists to prove. This is why
  `dual_hybrid` carries five records of each classification rather than one.
- **Keep `poison_gzip`'s SNP list small.** The failure is silent because a ~1 KB buffered write still
  fits the pipe buffer before the dead child is noticed. A larger list would behave differently.
- **Prefer error fixtures whose `die` carries no `$!`.** `strerror` text is not guaranteed identical
  across platforms.
- For `reverse_strand`, the reference base must be the *complement* of the file's ref allele, since
  `read_snps` complements before `create_modified_chromosome` compares.

## Normalisation

Beyond the shared masking: `*.txt.gz` is decompressed and recorded under a `.txt` suffix, following the
`.bam` → `.sam.txt` precedent, and `manifest` still pins the archive's existence separately from its
content.

**Six output files used to be written straight out of Perl hashes, so their line order differed on
every run** (#104). The source sorts them now, but the runner still sorts before comparison, so `--impl`
stays usable against an implementation with its own ideas about iteration order. The dual-hybrid
annotation additionally has its ID column's digits masked:

| Output | Cause |
|---|---|
| `all_SNPs_<strain>_<build>.txt.gz` | `keys %all_SNPs`, once per strain |
| `<strain>_specific_SNPs.<build>.txt` | `keys %snps` |
| `<strain2>_specific_SNPs.<build>.txt` | read order of the strain-2 archive |
| `<strain>_<strain2>_SNPs_in_common.<build>.txt` | read order of the strain-2 archive |
| `all_<strain2>_SNPs_<strain>_reference.based_on_<build>.txt` | line order **and** field 0 |
| the strain list on stdout | `keys %strains`, five call sites |

Sorting loses nothing: line order is read by no consumer, and `sort` keeps duplicates so a spurious
duplicate line is still caught. Field 0 is masked only in its digits — the `<strain>_` prefix stays, because
it is the one signal that a line came from the new-reference branch and it is a documented column.

**The table is hardcoded in the runner, not declared per fixture.** The nondeterminism is a property of
the tool, so an author who forgot a marker would get a fixture that fails one run in three, which is
worse than no fixture.

Perl stamps a line number onto the four `die`s and `warn`s that lack a trailing newline, so
`<impl> line NNN` is masked while `<IN> line N` — real behaviour — survives. Without it every fixture
pinning a diagnostic would break on any edit above the reported line.

Each fixture's own staged `genome/` and VCF are recorded into `expected/` **deliberately**: the tool
`chdir`s into the reference folder, so a regression that wrote there would otherwise be invisible.

## Fixtures worth knowing about

Four were written to pin behaviour that was wrong and now assert the fix, since #102, #103, #105 and
#106 are resolved:

- **`chrom_name_mismatch`** — an Ensembl-style VCF against a UCSC-style reference used to exit 0 and
  write `chrchr1.N-masked.fa` containing zero Ns. The per-chromosome check could not catch it, because
  the loop iterates the reference's chromosomes and skips any the VCF does not mention. Now checked once
  after both name sets are known, printing both lists. Partial overlap is still legal, which
  `multi_chrom` asserts.
- **`poison_gzip`** — a failing `gzip` used to leave a zero-byte archive and report success. `open` on a
  pipe returns before the shell runs, so the `close` is the only place it can be noticed.
- **`missing_genome`** — a missing `--reference_genome` used to exit 0.
- **`skip_filtering_dual`** — `--skip_filtering` used to accept `--strain2` and discard it, because both
  the promotion to `--dual_hybrid` and its implied `--full_sequence` sat inside the block
  `--skip_filtering` skips. It now builds all six genomes from committed SNP files with no VCF at all,
  which is the only fixture reaching `read_snp_files` against pre-made archives.

**`empty_chromosome`** reaches the Ensembl-versus-UCSC abort from the other direction: a header-only
FastA entry is stored as the empty string, which is falsy, so the per-chromosome check fires for a
chromosome that *is* present — reporting it as not found and then listing it among the names it found.
That path is unchanged and still worth knowing about.

**`no_format_column`** renames the `FORMAT` column rather than deleting it. Deleting it leaves nine
columns, `detect_strains` skips indices up to 8, and the run dies earlier on an empty strain list, which
would make the fixture a duplicate of `strain_not_in_vcf`.

**`duplicate_chrom_name`** reaches the end-of-file duplicate check, whose message ends in a full stop.
The in-loop check one screen above ends in an exclamation mark and needs a *third* occurrence of the
name, because it runs before the previous chromosome has been stored. The two messages are otherwise
identical.

**`ref_mismatch`** is the only fixture that deliberately writes a SNP annotation disagreeing with its
reference, so it is the one exception to the derive-never-type rule above — and the only way to reach
either skip counter, because a derived annotation always matches. It pins that 5 SNPs total, 3 applied,
1 already present and 1 mismatched all appear, rather than leaving two positions unaccounted for (#112).

**`bam_write_fails`** stages a samtools that runs normally and then reports failure when writing
genome1. Exiting straight away would be a different test: that kills `tag2sort` mid-write and was already
reported as `killed by signal 13`. The silent case — samtools finishing its work and reporting a problem,
which only `close` can observe — is what #116 was about, and the fixture records exit status 0 against
the previous version.

**`samtools_path_spaced`** and **`spaced_paths`** cover #99, one per suite. The first stages a samtools
wrapper at `samtools dir/samtools` and drives the whole pipeline through it; its `poison_calls` file is
empty, which is the assertion that `--samtools_path` was honoured at every site rather than just the
first. The second puts a gzipped VCF under `vcf dir/` and gives the build a name containing a space, so
it lands in the gzipped SNP list's file name too. Both produce nothing at all against the previous
version.

The spaced paths are substituted into the argument list *after* the `args` file has been split on
whitespace, via `<SPACED_SAMTOOLS>`, `<SPACED_VCF>` and `<SPACED_BUILD>`. A space written directly into
`args` would just become an argument boundary and test nothing.

**`reverse_strand`** is the only fixture that can reach the strand-complementing branch: VCF filtering
always writes strand `1`, so a `-1` can only arrive through a hand-written `--skip_filtering` annotation.

## Proving the fixtures assert anything

`29 passed` alone only proves the fixtures are self-consistent. Two checks establish more, both run by
hand:

```sh
# the normalisation is doing work, not getting lucky
for seed in 1 4242; do PERL_HASH_SEED=$seed PERL_PERTURB_KEYS=2 test/run_genome_tests.pl; done

# the tool itself produces byte-identical output on identical input
test/bin/check_reproducible.pl
```

CI runs both. The seed loop proves the normalised output is seed-independent; the ≥5-record rule above
is what makes the raw output differ when it should, and without that the check is vacuous.

`check_reproducible.pl` exists because **the suite cannot catch a reproducibility regression**: it sorts
those six outputs before comparing, so it passes whether or not the source sorts them. That script
compares the raw output of two runs under different hash seeds instead, and fails naming the files that
differ. Run against `SNPsplit_genome_preparation` as of `dev` before #104 was fixed, it names exactly
those six.

Mutation testing is the other half — change one line, confirm a *named* set of fixtures fails, revert:

| Mutate | Should fail |
|---|---|
| `'N'` in `create_modified_chromosome` | 14 fixtures, but **not** `chrom_name_mismatch`, `empty_chromosome` or `multi_chrom`'s VCF-absent chromosome, which legitimately introduce zero Ns |
| `'N'` in `create_modified_chromosome_dual_hybrid` | `dual_hybrid` only — the first mutation does not reach it |
| either full-sequence substitution | `dual_hybrid`, `full_sequence`, `no_nmasking` / `dual_hybrid` |
| `$fi == 1` | 18 fixtures |
| the strand complement in `read_snps` | `reverse_strand` only |
| the new Ref/SNP allele, the unique-to-strain-1 test, or the strain-1 filter test | `dual_hybrid` only |

A mutation that fails nothing means the code it touched is asserted by no fixture.

## Checking a release against real data

The fixtures cannot imitate a real library: millions of reads, every SNP shape, both strands. That is the
only thing that could have contradicted #94, where two bisulfite branches were removed on a control-flow
argument no fixture can reach.

```sh
test/bin/compare_versions.pl --bam sample.bam --snps all_SNPs_STRAIN_GRCm39.txt.gz \
    --old 0.8.0 --new dev --args '--paired --bisulfite'
```

It stages all three scripts at each revision — `SNPsplit` calls `$RealBin/tag2sort` with no override, so
mixing versions across that boundary would compare something nobody runs — runs both over the same BAM,
and compares alignment records in full.

BAM headers are compared with `@PG` removed, because `@PG` records the argv and 0.9.0 writes BAM through
`samtools view -o` rather than a shell redirect. `version`, `date_run` and `command` are filtered from
the reports for the same reason the YAML masks them.

Exit status is 0 when nothing that matters changed, 1 when something did, and **2 when either run
failed** — that last case exists because a run that produces nothing would otherwise pass the comparison
with no files to compare, which is the vacuous pass the rest of this suite is built to avoid.

## Environment

`SNPSPLIT_NO_SLEEP=1` skips the tools' progress pauses. Both harnesses set it; it takes a single run from
roughly 2.2 s to 0.2 s and matters most for `--verbose`, which otherwise sleeps once per MD element. All
three tools honour it.

The expected output in this repository was recorded with samtools 1.21 and perl 5.34.1 on macOS. CI runs
`ubuntu-24.04`. The same expected output has also been checked against samtools 1.22.1 with perl 5.32.1
and samtools 1.23.1 with perl 5.38.2 on Linux, and all three agree exactly — the samtools code that
writes SAM headers is unchanged from 1.19.1 to 1.21, and the volatile `@PG` and `@HD` fields are masked
anyway. If a future samtools does change what it writes, the fix is to re-record with `--update` and
check the diff is only the header line you expected.
