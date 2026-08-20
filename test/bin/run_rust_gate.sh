#!/bin/sh
# Run one fixture suite against the Rust implementation and compare the result against
# test/rust_fixtures.txt.
#
# Exits non-zero if a listed fixture fails (a regression) OR if an unlisted fixture passes
# (a stale list). The second case is not pedantry: it is what keeps a long port honest,
# because it makes claiming your own work mandatory rather than optional.
#
# A second argument is a worker count, exported as SNPSPLIT_PARALLEL so the same fixtures run
# against the same expected output at a different worker count. Output is supposed to be
# identical whatever the worker count, so this is that claim checked by the fixtures rather
# than asserted in a comment.
#
# Usage: test/bin/run_rust_gate.sh alignment|genome [workers]
set -eu

suite=${1:?usage: run_rust_gate.sh alignment|genome [workers]}
workers=${2:-}
root=$(cd "$(dirname "$0")/../.." && pwd)
cd "$root"

case "$suite" in
    alignment)
        runner="test/run_tests.pl"
        impl="rust/target/impl/SNPsplit"
        prefix=""
        all=$(ls test/fixtures)
        ;;
    genome)
        runner="test/run_genome_tests.pl"
        impl="rust/target/impl/SNPsplit_genome_preparation"
        prefix="genome:"
        all=$(ls test/genome_fixtures)
        ;;
    *)
        echo "unknown suite '$suite'" >&2
        exit 2
        ;;
esac

[ -x "$impl" ] || { echo "no Rust build at $impl; run rust/build-impl.sh first" >&2; exit 2; }

if [ -n "$workers" ]; then
    SNPSPLIT_PARALLEL=$workers
    export SNPSPLIT_PARALLEL
    echo "running with SNPSPLIT_PARALLEL=$workers"
fi

expected=$(grep -v '^[[:space:]]*#' test/rust_fixtures.txt | grep -v '^[[:space:]]*$' || true)

status=0
for fixture in $all; do
    listed=no
    for e in $expected; do
        [ "$e" = "$prefix$fixture" ] && listed=yes
    done

    if "$runner" --impl "$impl" "$fixture" >/dev/null 2>&1; then
        passed=yes
    else
        passed=no
    fi

    if [ "$listed" = yes ] && [ "$passed" = no ]; then
        echo "REGRESSION  $prefix$fixture is listed as passing but failed"
        status=1
    elif [ "$listed" = no ] && [ "$passed" = yes ]; then
        echo "STALE LIST  $prefix$fixture passes but is not listed; add it to test/rust_fixtures.txt"
        status=1
    elif [ "$listed" = yes ]; then
        echo "ok          $prefix$fixture"
    else
        echo "pending     $prefix$fixture"
    fi
done

exit $status
