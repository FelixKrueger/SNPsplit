#!/bin/sh
# Build the release binary and expose it under the three classic tool names, so the Perl
# fixture runners can be pointed at it:
#
#   test/run_tests.pl        --impl rust/target/impl/SNPsplit
#   test/run_genome_tests.pl --impl rust/target/impl/SNPsplit_genome_preparation
#
# The runners copy the implementation into a scratch directory before running it, so these
# are hard links rather than symlinks: a copied symlink would be a copy of its target under
# a name that still dispatches correctly, but a broken link on a filesystem without symlink
# support would fail confusingly. Hard links behave the same everywhere.
set -eu

cd "$(dirname "$0")"
cargo build --release --locked

out=target/impl
mkdir -p "$out"
for name in SNPsplit tag2sort SNPsplit_genome_preparation; do
    rm -f "$out/$name"
    ln target/release/snpsplit "$out/$name"
done

echo "built:"
ls -l "$out"
