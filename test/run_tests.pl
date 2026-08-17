#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Find ();
use File::Spec ();
use Cwd qw(abs_path getcwd);

### Differential regression harness for SNPsplit + tag2sort.
###
### Each fixture under test/fixtures/<name>/ is run in a private scratch directory and every file the
### run produces is normalised and compared against the committed files in <name>/expected/.
###
### Point --impl at any implementation (Perl today, Rust later) to diff it against the same expected output.

my $test_root = abs_path(dirname(__FILE__));
my $repo_root = abs_path(File::Spec->catdir($test_root, '..'));

my ($impl, $samtools, $update, $keep, $verbose, $help);

GetOptions(
    'update'     => \$update,
    'impl=s'     => \$impl,
    'samtools=s' => \$samtools,
    'keep'       => \$keep,
    'v|verbose'  => \$verbose,
    'h|help'     => \$help,
) or usage(1);

usage(0) if $help;

$impl     = defined $impl     ? abs_path($impl)     : File::Spec->catfile($repo_root, 'SNPsplit');
$samtools = defined $samtools ? abs_path($samtools) : which('samtools');

die "SNPsplit implementation not found or not executable: $impl\n" unless -x $impl;
die "samtools not found; install it or pass --samtools PATH\n"     unless defined $samtools && -x $samtools;

my $tag2sort = File::Spec->catfile(dirname($impl), 'tag2sort');
die "tag2sort not found beside the implementation: $tag2sort\n" unless -x $tag2sort;

my $fixture_dir = File::Spec->catdir($test_root, 'fixtures');
my @fixtures    = @ARGV ? @ARGV : all_fixtures();
die "No fixtures found under $fixture_dir\n" unless @fixtures;

### Named up front so a failure in a log is attributable without reproducing the run: --verbose alone
### only ever shows the per-fixture staged copy, which is a fresh temporary path every time.
printf "Testing %s\nUsing   %s\n\n", $impl, $samtools;

my (@passed, @failed);

for my $name (@fixtures) {
    my $dir = File::Spec->catdir($fixture_dir, $name);
    unless (-d $dir) {
        warn "FAIL $name (no such fixture directory)\n";
        push @failed, $name;
        next;
    }
    my $ok = run_fixture($name, $dir);
    $ok ? push @passed, $name : push @failed, $name;
}

print "\n";
printf "%s: %d passed, %d failed (of %d)\n",
    ($update ? 'Recorded' : 'Result'), scalar @passed, scalar @failed, scalar @fixtures;
print "Failed: @failed\n" if @failed;

exit(@failed ? 1 : 0);

###############################################################################
### Fixture execution
###############################################################################

sub run_fixture {
    my ($name, $dir) = @_;

    my $scratch = tempdir("snpsplit_${name}_XXXXXX",
                          DIR     => ($ENV{TMPDIR} || '/tmp'),
                          CLEANUP => !$keep);
    print "  scratch: $scratch\n" if $verbose || $keep;

    ### $RealBin in SNPsplit resolves symlinks and hardcodes "$RealBin/tag2sort", so the pair must be
    ### staged together as real copies for --impl to be honoured.
    my $bin = File::Spec->catdir($scratch, 'bin');
    make_path($bin);
    my $staged_impl = File::Spec->catfile($bin, 'SNPsplit');
    copy($impl, $staged_impl)                                  or die "copy $impl: $!\n";
    my $staged_sorter = File::Spec->catfile($bin, 'tag2sort');
    copy($tag2sort, $staged_sorter)                            or die "copy $tag2sort: $!\n";

    ### Replacing the staged sorter with one that always fails tests the contract "SNPsplit aborts when
    ### the sorting stage fails" on its own terms. Relying on a real failure means the coverage
    ### disappears the moment that failure is fixed.
    if (-e File::Spec->catfile($dir, 'break_tag2sort')) {
        open my $fh, '>', $staged_sorter or die "write failing sorter: $!\n";
        print {$fh} "#!/bin/sh\necho 'tag2sort: deliberate failure' >&2\nexit 3\n";
        close $fh;
    }
    chmod 0755, $staged_impl, $staged_sorter;

    copy(File::Spec->catfile($dir, 'snps.txt'), File::Spec->catfile($scratch, 'snps.txt'))
        or die "$name: fixture has no readable snps.txt: $!\n";

    my @args     = read_args($dir);
    my $feed_sam = -e File::Spec->catfile($dir, 'feed_sam');
    my $poison   = -e File::Spec->catfile($dir, 'poison_shim');

    my @inputs = prepare_inputs($name, $dir, $scratch, $feed_sam);
    die "$name: no input*.sam found\n" unless @inputs;

    ### The poison shim proves --samtools_path is honoured everywhere: any bare `samtools` call lands
    ### on the shim, which records itself and fails. PATH subtraction cannot be used instead, because
    ### PATH must still resolve perl (both scripts are #!/usr/bin/env perl) and on Debian-family CI
    ### samtools and perl share /usr/bin.
    my $run_path = $ENV{PATH};
    my $poison_log;
    if ($poison) {
        my $shim_dir = File::Spec->catdir($scratch, 'shim');
        make_path($shim_dir);
        $poison_log = File::Spec->catfile($scratch, 'poison_calls');
        my $shim = File::Spec->catfile($shim_dir, 'samtools');
        open my $fh, '>', $shim or die "write shim: $!\n";
        print {$fh} <<"SHIM";
#!/bin/sh
echo "bare samtools invoked: \$*" >> "$poison_log"
echo "samtools: shim reached - --samtools_path was ignored" >&2
exit 127
SHIM
        close $fh;
        chmod 0755, $shim;
        $run_path = "$shim_dir:$run_path";
        open my $touch, '>', $poison_log or die "create poison log: $!\n";
        close $touch;
    }

    s/<REAL_SAMTOOLS>/$samtools/g for @args;

    ### Relative filenames are mandatory: the YAML records infile and SNP_annotation verbatim, so an
    ### absolute path would bake the random scratch component into the recorded output.
    my @cmd = ($staged_impl, '--SNP_file', 'snps.txt', @args, @inputs);
    my $cmdline = join ' ', map { shell_quote($_) } @cmd;
    print "  run: $cmdline\n" if $verbose;

    my $prev = getcwd();
    chdir $scratch or die "chdir $scratch: $!\n";
    local $ENV{PATH}             = $run_path;
    local $ENV{SNPSPLIT_NO_SLEEP} = 1;
    my $status = system("$cmdline > stdout.log 2> stderr.log");
    my $exit   = $status == -1 ? -1 : $status >> 8;
    chdir $prev or die "chdir back: $!\n";

    my %actual = normalise($name, $dir, $scratch, $exit);

    return $update
        ? write_expected($name, $dir, \%actual, \@inputs)
        : compare_expected($name, $dir, \%actual);
}

sub prepare_inputs {
    my ($name, $dir, $scratch, $feed_sam) = @_;

    my @src = sort glob(File::Spec->catfile($dir, 'input*.sam'));
    my @inputs;
    my $n = 0;
    for my $src (@src) {
        check_md_tags($name, $dir, $src);
        $n++;
        my $stem = $n == 1 ? $name : "${name}_$n";
        if ($feed_sam) {
            my $dst = "$stem.sam";
            copy($src, File::Spec->catfile($scratch, $dst)) or die "copy $src: $!\n";
            push @inputs, $dst;
        }
        else {
            my $dst = File::Spec->catfile($scratch, "$stem.bam");
            my $st  = system(join ' ',
                shell_quote($samtools), 'view', '-bS', shell_quote($src), '>', shell_quote($dst));
            ### An unchecked conversion lets --update record plausible all-empty output from a
            ### mis-authored fixture.
            die "$name: samtools could not convert $src (exit @{[ $st >> 8 ]}) - check CIGAR/SEQ lengths\n"
                if $st != 0;
            push @inputs, "$stem.bam";
        }
    }
    return @inputs;
}

### The expected output is only as trustworthy as the MD tags, and nothing in samtools validates one:
### a tag that is self-consistent but names the wrong coordinate would pin the wrong behaviour and look
### exactly like a correct fixture. Re-deriving from the committed reference means a hand-edited
### input.sam cannot get that far.
sub check_md_tags {
    my ($name, $dir, $src) = @_;

    my $fa = File::Spec->catfile($dir, 'chr1.fa');
    return unless -e $fa;

    my %committed;
    for my $line (split /^/, slurp($src)) {
        next if $line =~ /^\@/;
        my @f = split /\t/, $line;
        my ($md) = map { /^MD:Z:(\S+)/ ? $1 : () } @f[11 .. $#f];
        $committed{ $f[0] . "\t" . $f[3] } = $md if defined $md;
    }
    return unless %committed;

    ### calmd wants a faidx-able reference and writes its own index beside it.
    my $q_fa  = shell_quote($fa);
    my $q_src = shell_quote($src);
    system("$samtools faidx $q_fa 2>/dev/null") == 0
        or die "$name: could not index $fa\n";
    my $rederived = `$samtools calmd $q_src $q_fa 2>/dev/null`;
    unlink "$fa.fai";
    die "$name: calmd could not re-derive MD tags from chr1.fa\n" if $? != 0 || !length $rederived;

    for my $line (split /^/, $rederived) {
        next if $line =~ /^\@/;
        my @f = split /\t/, $line;
        my ($md) = map { /^MD:Z:(\S+)/ ? $1 : () } @f[11 .. $#f];
        next unless defined $md;
        my $key = $f[0] . "\t" . $f[3];
        next unless exists $committed{$key};       # reads with MD deliberately stripped
        next if $committed{$key} eq $md;
        die "$name: MD tag in " . basename($src) . " does not match chr1.fa\n"
          . "  read $f[0] at position $f[3]\n"
          . "  committed:  MD:Z:$committed{$key}\n"
          . "  reference:  MD:Z:$md\n"
          . "  Re-author with test/bin/make_fixtures.pl rather than editing input.sam by hand.\n";
    }
}

###############################################################################
### Normalisation
###############################################################################

sub normalise {
    my ($name, $dir, $scratch, $exit) = @_;

    my @excludes = read_lines(File::Spec->catfile($dir, 'exclude'));
    my @present  = collect_files($scratch);

    my %out;
    ### Only whether the run aborted, not the number. Perl's die exits with $!, which reflects the last
    ### failed syscall and therefore varies with whatever ran before; the die message itself is pinned
    ### in run.log, so the number carried no information and made the fixture intermittent.
    $out{'exit_status'} = $exit == 0 ? "0\n" : "nonzero\n";
    $out{'manifest'}    = join('', map { "$_\n" } @present);

    for my $rel (@present) {
        next if $rel eq 'stdout.log' || $rel eq 'stderr.log';
        next if excluded($rel, \@excludes);

        my $abs = File::Spec->catfile($scratch, $rel);

        if ($rel =~ /\.bam$/) {
            $out{"$rel.sam.txt"} = bam_to_text($abs);
        }
        elsif ($rel =~ /\.sam$/) {
            $out{"$rel.txt"} = mask_sam_text(slurp($abs));
        }
        elsif ($rel =~ /\.yaml$/) {
            $out{$rel} = mask_yaml(slurp($abs));
        }
        elsif ($rel eq 'poison_calls') {
            ### Concurrent samtools processes append here, so the order of the lines is a scheduling
            ### artifact: sam2bam_fails records two calls and failed roughly one run in four. The set of
            ### calls is what the fixture asserts, so sorting loses nothing.
            my @calls = sort split /^/, mask_paths(slurp($abs), $scratch);
            $out{$rel} = join '', @calls;
        }
        else {
            $out{$rel} = mask_paths(slurp($abs), $scratch);
        }
    }

    ### A third of SNPsplit's user-visible output is stderr-only - the hard-clip counter at
    ### SNPsplit:484 among it - so stderr is recorded and compared, not just a debugging aid.
    $out{'run.log'} = mask_paths(slurp(File::Spec->catfile($scratch, 'stderr.log')), $scratch);
    my $stdout = slurp(File::Spec->catfile($scratch, 'stdout.log'));
    $out{'stdout.log'} = mask_paths($stdout, $scratch) if length $stdout;

    return %out;
}

sub bam_to_text {
    my ($abs) = @_;

    ### The suite must survive judging deliberately-broken output, so an unreadable BAM becomes a
    ### comparable sentinel rather than aborting the run.
    return "<unreadable>\n" if -z $abs;

    my $q      = shell_quote($abs);
    my $text   = `@{[ shell_quote($samtools) ]} view -h $q 2>/dev/null`;
    my $failed = $? != 0;
    return "<unreadable>\n" if $failed || !length $text;

    return mask_sam_text($text);
}

sub mask_sam_text {
    my ($text) = @_;
    my @out;
    for my $line (split /^/, $text) {
        ### @PG is a control input, not decoration: check_for_bs reads it to decide library type and
        ### bisulfite mode, so the IDs and PP: chain stay and only the volatile fields are masked.
        if ($line =~ /^\@PG/) {
            $line =~ s/\tVN:[^\t\n]*/\tVN:<masked>/g;
            $line =~ s/\tCL:[^\t\n]*/\tCL:<masked>/g;
        }
        ### SS: is a recent sub-sort tag whose presence varies by samtools version; SO: is behaviour.
        elsif ($line =~ /^\@HD/) {
            $line =~ s/\tSS:[^\t\n]*/\tSS:<masked>/g;
        }
        push @out, $line;
    }
    return join '', @out;
}

sub mask_yaml {
    my ($text) = @_;
    for my $key (qw(version date_run command)) {
        $text =~ s/^(\s*\Q$key\E:).*$/$1 <masked>/mg;
    }
    return $text;
}

sub mask_paths {
    my ($text, $scratch) = @_;
    my $real = defined $scratch ? abs_path($scratch) : undef;
    for my $p (grep { defined && length } ($real, $scratch)) {
        $text =~ s/\Q$p\E/<scratch>/g;
    }
    ### Both tools resolve samtools with `which`, so they report the PATH entry while the runner holds
    ### the symlink-resolved path. Masking the reported value keeps the install location out of the
    ### recorded output whatever it is; matching on $samtools alone leaves it baked in.
    $text =~ s/^(Samtools path:\s*)\S+/$1<samtools>/mg;
    $text =~ s/\Q$samtools\E/<samtools>/g;
    return $text;
}

sub collect_files {
    my ($scratch) = @_;
    my @files;
    File::Find::find(
        {
            no_chdir => 1,
            wanted   => sub {
                my $path = $File::Find::name;
                return unless -f $path;
                my $rel = File::Spec->abs2rel($path, $scratch);
                ### bin/ and shim/ are harness scaffolding, not results.
                return if $rel =~ m{^(?:bin|shim)/};
                push @files, $rel;
            },
        },
        $scratch
    );
    return sort @files;
}

###############################################################################
### Comparison
###############################################################################

sub write_expected {
    my ($name, $dir, $actual, $inputs) = @_;

    ### Validation 2 belongs inside --update: a fixture whose outputs are all empty is almost always a
    ### mis-authored fixture, and recording it would pin the mistake. The fixture's own input
    ### is excluded from the check, or its records satisfy it single-handed and every empty output walks.
    ### It runs before expected/ is touched, so a refusal leaves the previous files intact.
    my %is_input    = map { ("$_.sam.txt" => 1, "$_.txt" => 1) } @$inputs;
    my @outputs     = grep { /\.sam\.txt$/ && !$is_input{$_} } keys %$actual;
    my $allow_empty = -e File::Spec->catfile($dir, 'allow_empty');

    my @unreadable = sort grep { $actual->{$_} =~ /^<unreadable>/ } @outputs;
    if (@unreadable && !$allow_empty) {
        warn "FAIL $name: alignment output could not be read (@unreadable); refusing to record it.\n"
           . "      If that is genuinely intended, add an 'allow_empty' marker to the fixture.\n";
        return 0;
    }

    my $any_records = 0;
    for my $k (@outputs) {
        $any_records = 1 if grep { !/^\@/ && /\S/ } split /^/, $actual->{$k};
    }
    if (@outputs && !$any_records && !$allow_empty) {
        warn "FAIL $name: every alignment output is empty; refusing to record it.\n"
           . "      If that is genuinely intended, add an 'allow_empty' marker to the fixture.\n";
        return 0;
    }

    my $exp = File::Spec->catdir($dir, 'expected');
    if (-d $exp) {
        for my $f (glob File::Spec->catfile($exp, '*')) { unlink $f or die "unlink $f: $!\n" }
    }
    else { make_path($exp) }

    for my $rel (sort keys %$actual) {
        my $dst = File::Spec->catfile($exp, flatten($rel));
        open my $fh, '>', $dst or die "write $dst: $!\n";
        print {$fh} $actual->{$rel};
        close $fh;
    }
    printf "UPDATED %-20s %d expected files\n", $name, scalar keys %$actual;
    return 1;
}

sub compare_expected {
    my ($name, $dir, $actual) = @_;

    my $exp = File::Spec->catdir($dir, 'expected');
    unless (-d $exp) {
        warn "FAIL $name: no expected/ directory - run with --update first\n";
        return 0;
    }

    my %expected;
    for my $f (glob File::Spec->catfile($exp, '*')) {
        $expected{ basename($f) } = slurp($f);
    }
    my %got = map { flatten($_) => $actual->{$_} } keys %$actual;

    my @diffs;
    for my $k (sort keys %expected) {
        if (!exists $got{$k})            { push @diffs, "  missing output: $k" }
        elsif ($got{$k} ne $expected{$k})  { push @diffs, "  differs: $k", indent(diff_text($expected{$k}, $got{$k})) }
    }
    for my $k (sort keys %got) {
        push @diffs, "  unexpected output: $k" unless exists $expected{$k};
    }

    if (@diffs) {
        print "FAIL $name\n", join("\n", @diffs), "\n";
        return 0;
    }
    printf "PASS %s\n", $name;
    return 1;
}

### Minimal line diff; keeps the failure readable without a diff(1) dependency or TAP indirection.
sub diff_text {
    my ($want, $got) = @_;
    my @w = split /^/, $want;
    my @g = split /^/, $got;
    my @out;
    my $max = @w > @g ? @w : @g;
    for my $i (0 .. $max - 1) {
        my $a = defined $w[$i] ? $w[$i] : '';
        my $b = defined $g[$i] ? $g[$i] : '';
        next if $a eq $b;
        chomp(my $ca = $a);
        chomp(my $cb = $b);
        push @out, sprintf("line %d:\n  -expected: %s\n  +actual:   %s", $i + 1, $ca, $cb);
        last if @out >= 10;
    }
    push @out, '(further differences suppressed)' if @out >= 10;
    return join "\n", @out;
}

###############################################################################
### Helpers
###############################################################################

sub all_fixtures {
    return sort map { basename($_) } grep { -d } glob File::Spec->catfile($fixture_dir, '*');
}

sub read_args {
    my ($dir) = @_;
    my @lines = read_lines(File::Spec->catfile($dir, 'args'));
    return map { split ' ' } @lines;
}

sub read_lines {
    my ($file) = @_;
    return () unless -e $file;
    open my $fh, '<', $file or die "read $file: $!\n";
    my @lines = grep { /\S/ && !/^\s*#/ } map { chomp; $_ } <$fh>;
    close $fh;
    return @lines;
}

sub slurp {
    my ($file) = @_;
    return '' unless -e $file;
    open my $fh, '<', $file or return '';
    local $/;
    my $c = <$fh>;
    close $fh;
    return defined $c ? $c : '';
}

sub excluded {
    my ($rel, $patterns) = @_;
    for my $p (@$patterns) {
        my $re = quotemeta $p;
        $re =~ s/\\\*/.*/g;
        return 1 if $rel =~ /^$re$/;
    }
    return 0;
}

### expected/ is one flat directory, so nested output paths (e.g. --output_dir) are encoded.
sub flatten {
    my ($rel) = @_;
    $rel =~ s{/}{__}g;
    return $rel;
}

sub indent { my ($t) = @_; $t =~ s/^/    /mg; return $t }

sub shell_quote {
    my ($s) = @_;
    return $s if $s =~ m{^[A-Za-z0-9_./:=-]+$};
    $s =~ s/'/'\\''/g;
    return "'$s'";
}

sub which {
    my ($prog) = @_;
    for my $d (split /:/, ($ENV{PATH} || '')) {
        my $p = File::Spec->catfile($d, $prog);
        return abs_path($p) if -x $p && !-d $p;
    }
    return undef;
}

sub usage {
    my ($code) = @_;
    print <<'USAGE';
Differential regression harness for SNPsplit + tag2sort.

  test/run_tests.pl [OPTIONS] [FIXTURE ...]

  --update          Regenerate expected/ for the selected fixtures instead of diffing.
  --impl PATH       SNPsplit executable under test. Default: the repo's ./SNPsplit.
                    tag2sort is taken from the same directory and staged alongside it.
  --samtools PATH   samtools for fixture prep and normalisation. Default: first on PATH.
  --keep            Leave scratch directories in place for inspection.
  -v, --verbose     Echo each command as it runs.
  -h, --help        This message.

  With no FIXTURE arguments, every directory under test/fixtures/ runs.
  Exit status is 0 when all selected fixtures match, 1 otherwise.

  Never run --update without reading the resulting diff: it will happily pin a bug.
USAGE
    exit $code;
}
