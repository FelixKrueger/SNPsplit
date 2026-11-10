package FixedClock;
# Loaded through PERL5OPT by test/run_tests.pl: localtime() is noon UTC on 6 Nov 2026, a one-digit day
use strict;
use warnings;
BEGIN { *CORE::GLOBAL::localtime = sub (;$) { CORE::localtime(1793966400) } }
1;
