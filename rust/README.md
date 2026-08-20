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
| `genome_prep` | not started | 0 / 31 |
| `sort` (tag2sort) | not started | 0 / 26 |
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

None yet. Anything that lands here stays here: a deviation that is not written down is a
bug the next person has to rediscover.
