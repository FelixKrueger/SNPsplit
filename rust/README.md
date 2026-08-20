# SNPsplit Rust rewrite

Port of SNPsplit from Perl to Rust, proposed in
[#89](https://github.com/FelixKrueger/SNPsplit/issues/89). The design is in
[DESIGN.md](DESIGN.md).

This file is the status journal. It says what is ported and at what fidelity, and it is
updated in the same PR as the work it describes.

## Building

```sh
./rust/build-impl.sh
```

Builds one binary and exposes it under the three classic names in `rust/target/impl/`,
which is what the fixture runners take:

```sh
test/run_tests.pl        --impl rust/target/impl/SNPsplit
test/run_genome_tests.pl --impl rust/target/impl/SNPsplit_genome_preparation
```

## Status

| Module | State | Fixtures passing |
|---|---|---|
| multicall dispatch, version banners | done | n/a |
| `io` (BAM/SAM read, BAM write, name sort) | done | n/a |
| `genome_prep` | done | 31 / 31 |
| `sort` (tag2sort) | done | 25 / 26 driven by the Perl tagger |
| `tag` (SNPsplit) | not started | 0 / 26 |

The two fixture counts are the two suites, not two halves of one: the alignment suite runs
`SNPsplit` and `tag2sort` together, so its 26 fixtures only pass once both are ported.

Which fixtures the Rust build is expected to pass is not a claim in this table but a
checked file, `test/rust_fixtures.txt`, enforced by `test/bin/run_rust_gate.sh` in both
directions.

## Verified equivalences

- All three version banners are byte-identical to Perl v0.9.0, checked by `diff` against the
  Perl scripts and pinned by a test that fails if `rust/VERSION` moves without the banners
  being regenerated.
- Name sort: `io::sort_by_name` reproduces `samtools sort -n` ordering, cross-checked over
  5000 shuffled names with a 500-record memory budget (so the spill and merge paths are the
  ones under test) against samtools 1.24. The comparator is a port of samtools' `strnum_cmp`,
  not a byte-wise string compare: digit runs compare by value, and equal values with
  different zero padding are ordered by the padding rather than treated as equal, so
  `read007` precedes `read7`.
- Our BAM output is readable by samtools, asserted directly rather than by comparing BGZF
  bytes. What the fixtures diff is samtools-rendered SAM text, so that is the contract.
- Worker-count invariance: the sorted output is identical at 1, 2 and 8 workers, and
  identical across memory budgets of 64, 333 and 100000 records. Both are tests, not claims.

## Measured

Name sort of 2,000,000 shuffled 50bp reads (18 MB BAM), 400k records in memory, on a
16-core machine:

| workers | elapsed |
|---|---|
| 1 | 4.29s |
| 2 | 3.05s |
| 4 | 2.45s |
| 8 | 2.17s |
| 16 | 2.06s |

`samtools sort -n` on the same input: 2.59s at one thread, 0.93s at eight. We are still
slower, and the remaining gap is the per-record write path, not the sort. Recorded here
rather than omitted: the point of the port is dropping the samtools dependency, not beating
it.

Scaling flattens after four workers because the merge is sequential by design. Before the
raw-record path landed, the single-worker time was 9.28s; decoding every record into an
owned `RecordBuf` only to re-encode it was two thirds of the cost.

## Known deviations from Perl v0.9.0

Anything that lands here stays here: a deviation that is not written down is a bug the next
person has to rediscover.

- **The all-SNP list is compressed through `gzip -c` when there is one to run**, and in
  process when there is not. Compressing in process throughout would be the obvious choice
  and produces bytes that decompress to the same text. The reason not to is that a *failing*
  gzip is observable behaviour: SNPsplit reads this file, and a truncated one means it loads
  no SNPs at all, so a gzip failure is fatal and `poison_gzip` pins that. With no subprocess
  there is nothing for the fixture's shim to poison. The fallback keeps the binary working
  where the Perl would not run at all. The two paths differ in which failures are possible,
  never in content. Reading a gzipped VCF is in process unconditionally: no fixture depends on
  that being a subprocess, and not forking is strictly better.
- **Perl `die` suffixes are emulated.** Four genome fixtures diff messages that Perl stamps
  with ` at <script> line <n>, <IN> line <m>.` because the `die` string has no trailing
  newline. The script line number is masked by the runner and carries no information; the
  input line number is not masked and is genuinely useful. The Rust build prints the same
  shape, using the Perl line numbers as provenance. The real fix is a trailing newline in the
  three `die` calls, which is raised upstream separately.

### The `@PG` chain

Output files record who wrote them. The Perl pipeline records only the `samtools` invocations
it shells out to, so the tool that actually did the work never appears in its own output. A
self-contained build has no excuse for that, so `tag2sort` writes its own record:

```
@PG	ID:SNPsplit	PN:SNPsplit	VN:0.9.0	CL:<the command line>	PP:<the previous record>
```

This is the one place the port writes something the Perl does not, and it is deliberate: the
alternative that would have matched byte for byte was emitting `PN:samtools` records for
invocations that never happened, which is false provenance in a data file.

The fixture runner now drops the provenance record an implementation's own I/O layer adds,
on either side, for the same reason it already masked those records' `VN:` and `CL:` fields.
Every other `@PG` record, including the aligner's that `check_for_bs` reads, is compared as
before.

### `bam_write_fails` is not listed

It poisons `samtools` so that writing `genome1.bam` fails, and asserts the run aborts rather
than reporting success over a truncated file. With no samtools to poison the write succeeds,
and the abort message it pins names samtools explicitly:

```
samtools failed while writing:
  bam_write_fails.genome1.bam
These output files are incomplete
```

Printing that from a build that never ran samtools is the same false claim as the `@PG`
records. Re-expressing the fixture needs the Perl message to drop the tool name so both
implementations can share it, which is raised upstream.

## Perl behaviour reproduced deliberately

These look like port bugs in a diff and are not. Each is raised upstream on its own.

- An **empty chromosome** is reported as a chromosome-name mismatch. `create_modified_chromosome`
  guards on `unless ($chromosomes{$chr})`, and an empty string is false in Perl, so a
  chromosome that exists but carries no sequence takes the "not found in the reference
  genome" path.
- A **chromosome present in the VCF body but never declared in a `##contig` header** aborts
  the run with `Can't use an undefined value as a symbol reference`: there is no filehandle
  for it, and the Perl reaches `print {$fhs{$chr}}` with an undefined value.
- **`Clearing SNP array...`** is printed for every chromosome with no SNP file. It announces
  assigning an empty list to a list that is already empty.
- The **dual hybrid report drops separators** that its stderr copy carries, differently on two
  adjacent lines: the N-masked line says `strainstrain 2 [AB]` and the full-sequence line says
  `strainstrain 2 [A/B]`, where stderr says `strain/strain 2 [A/B]` for both.
- **Duplicate SNP positions are counted before being skipped**, so the per-chromosome total
  does not reconcile with the applied and skipped counts.
