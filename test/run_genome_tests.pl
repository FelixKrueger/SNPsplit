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

### Differential regression harness for SNPsplit_genome_preparation.
###
### Each fixture under test/genome_fixtures/<name>/ is run in a private scratch directory and every file
### the run produces is normalised and compared against the committed files in <name>/expected/.
###
### Sibling of run_tests.pl, which drives SNPsplit and tag2sort. The tool-agnostic subs below are copied
### byte-identical from it so the two can be de-duplicated into a shared module mechanically; the four
### that could not be copied (mask_paths, write_expected, all_fixtures, usage) are noted at their
### definitions.
###
### Unlike run_tests.pl this driver injects no arguments: each fixture's args file is the complete
### command line. This tool takes three independent mandatory-ish inputs rather than one, so injection
### would leave a reader guessing which options the runner owns.

my $test_root = abs_path(dirname(__FILE__));
my $repo_root = abs_path(File::Spec->catdir($test_root, '..'));

my ($impl, $update, $keep, $verbose, $help);

GetOptions(
    'update'    => \$update,
    'impl=s'    => \$impl,
    'keep'      => \$keep,
    'v|verbose' => \$verbose,
    'h|help'    => \$help,
) or usage(1);

usage(0) if $help;

$impl = defined $impl ? abs_path($impl)
                      : File::Spec->catfile($repo_root, 'SNPsplit_genome_preparation');

die "Implementation not found or not executable: $impl\n" unless -x $impl;

### gzip and gunzip are unconditional runtime requirements of the tool, not just of the harness: it
### writes the all-SNP file through `| gzip -c` and reads gzipped VCFs through `gunzip -c |`. samtools
### is deliberately absent - nothing in this tool touches an alignment.
for my $prog (qw(gzip gunzip)) {
    die "$prog not found on PATH; the tool requires it\n" unless defined which($prog);
}

###############################################################################
### Normalisation table
###############################################################################

### Six output files carry hash iteration order, so byte comparison would make their fixtures fail
### intermittently. The table is hardcoded rather than declared per fixture on purpose: the
### nondeterminism is a property of the tool, and a fixture author who forgot the marker would get an
### intermittent fixture, which is worse than no fixture.
###
### Sorting is lossless here. Line order carries no information any consumer reads - SNPsplit:1431 and
### read_new_snp_annotation:317 both index by field - and sort preserves duplicates, so a spurious
### duplicate line is still caught.
###
### Declared before the fixture loop, not beside the other normalisation subs: a file-scoped assignment
### below the loop would still be empty when the loop ran.
my @NORMALISE = (
    ### keys %all_SNPs at :1196, once per strain
    { pattern => qr/^all_SNPs_.*\.txt$/,                           sort_lines => 1 },
    ### strain 1 from keys %snps at :510; strain 2 from the read order of the strain-2 archive at
    ### :497 - two different causes, one fix
    { pattern => qr/_specific_SNPs\..*\.txt$/,                     sort_lines => 1 },
    ### read order of the strain-2 archive at :449
    { pattern => qr/_SNPs_in_common\..*\.txt$/,                    sort_lines => 1 },
    ### Line order plus field 0: :465 counts over a hash-ordered file, :499 passes a line through
    ### verbatim, and :538 numbers with ${strain}_${unique_ref} in %snps order. The three collide -
    ### two different lines can carry ID 3 in one run - so the digits are masked and the prefix kept.
    { pattern => qr/^all_.*_SNPs_.*_reference\.based_on_.*\.txt$/, sort_lines => 1, mask_id => 1 },
);

my $fixture_dir = File::Spec->catdir($test_root, 'genome_fixtures');
my @fixtures    = @ARGV ? @ARGV : all_fixtures($fixture_dir);
die "No fixtures found under $fixture_dir\n" unless @fixtures;

printf "Testing %s\n\n", $impl;

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

    my $scratch = tempdir("gp_${name}_XXXXXX",
                          DIR     => ($ENV{TMPDIR} || '/tmp'),
                          CLEANUP => !$keep);
    print "  scratch: $scratch\n" if $verbose || $keep;

    my $bin = File::Spec->catdir($scratch, 'bin');
    make_path($bin);
    my $staged_impl = File::Spec->catfile($bin, basename($impl));
    copy($impl, $staged_impl) or die "copy $impl: $!\n";
    chmod 0755, $staged_impl;

    stage_inputs($name, $dir, $scratch);

    ### The poison shim proves nothing about --samtools_path here; it proves that a failed gzip is
    ### silent. PATH subtraction cannot be used instead, because PATH must still resolve perl.
    my $run_path = $ENV{PATH};
    if (-e File::Spec->catfile($dir, 'poison_gzip')) {
        my $shim_dir = File::Spec->catdir($scratch, 'shim');
        make_path($shim_dir);
        my $shim = File::Spec->catfile($shim_dir, 'gzip');
        open my $fh, '>', $shim or die "write shim: $!\n";
        print {$fh} "#!/bin/sh\necho 'gzip: deliberate failure' >&2\nexit 9\n";
        close $fh;
        chmod 0755, $shim;
        $run_path = "$shim_dir:$run_path";
    }

    my @args = read_args($dir);
    die "$name: fixture has no args file\n" unless @args;

    ### Paths containing a space, substituted after the args file has been split on whitespace so the
    ### space cannot become an argument boundary here. A VCF under a directory with a space in the name
    ### is the realistic case - the tool reads it through gunzip - and a build name with a space lands
    ### in every output file name, including the gzipped SNP list.
    if (grep { /<SPACED_VCF>/ } @args) {
        my $spaced = File::Spec->catdir($scratch, 'vcf dir');
        make_path($spaced);
        my @staged = grep { -f } ( glob(File::Spec->catfile($scratch, '*.vcf')),
                                   glob(File::Spec->catfile($scratch, '*.vcf.gz')) );
        die "$name: <SPACED_VCF> needs exactly one staged VCF, found @{[ scalar @staged ]}\n"
            unless @staged == 1;
        my $moved = File::Spec->catfile($spaced, basename($staged[0]));
        rename $staged[0], $moved or die "move VCF into '$spaced': $!\n";
        my $rel = File::Spec->abs2rel($moved, $scratch);
        s/<SPACED_VCF>/$rel/g for @args;
    }
    s/<SPACED_BUILD>/test build/g for @args;

    my @cmd     = ($staged_impl, @args);
    my $cmdline = join ' ', map { shell_quote($_) } @cmd;
    print "  run: $cmdline\n" if $verbose;

    my $prev = getcwd();
    chdir $scratch or die "chdir $scratch: $!\n";
    local $ENV{PATH}              = $run_path;
    local $ENV{SNPSPLIT_NO_SLEEP} = 1;
    my $status = system("$cmdline > stdout.log 2> stderr.log");
    my $exit   = $status == -1 ? -1 : $status >> 8;
    chdir $prev or die "chdir back: $!\n";

    my %actual = normalise($name, $dir, $scratch, $exit, basename($impl));

    return $update
        ? write_expected($name, $dir, \%actual, $exit)
        : compare_expected($name, $dir, \%actual);
}

### The reference is a directory rather than a file, which is the main departure from run_tests.pl's
### prepare_inputs. Staged inputs are collected into expected/ deliberately: the tool chdirs into the
### reference folder at :1283, so a regression that wrote there would otherwise be invisible.
sub stage_inputs {
    my ($name, $dir, $scratch) = @_;

    ### Optional: the fixtures that abort inside process_commandline never reach the genome, and
    ### shipping one would only add noise to their manifest.
    my $genome = File::Spec->catdir($dir, 'genome');
    copy_tree($genome, File::Spec->catdir($scratch, 'genome')) if -d $genome;

    ### Empty files created before the run, to reach the two branches that report finding their own
    ### output already in place.
    for my $rel (read_lines(File::Spec->catfile($dir, 'pre_existing'))) {
        my $abs = File::Spec->catfile($scratch, $rel);
        make_path(dirname($abs));
        open my $fh, '>', $abs or die "create $abs: $!\n";
        close $fh;
    }

    ### snps_in/ holds a pre-made SNPs_<strain>/ tree for --skip_filtering fixtures.
    my $snps_in = File::Spec->catdir($dir, 'snps_in');
    if (-d $snps_in) {
        opendir my $dh, $snps_in or die "opendir $snps_in: $!\n";
        for my $entry (grep { $_ ne '.' && $_ ne '..' } readdir $dh) {
            my $src = File::Spec->catfile($snps_in, $entry);
            if (-d $src) {
                copy_tree($src, File::Spec->catdir($scratch, $entry));
            }
            ### The all-SNP archives are stored as plain text and compressed on the way in, so no
            ### binary lands in the repository and regenerating a fixture cannot churn the diff with a
            ### different gzip version's output.
            elsif ($entry =~ /\.gz$/) {
                my $dst = File::Spec->catfile($scratch, $entry);
                my $st  = system(join ' ', 'gzip', '-c', shell_quote($src), '>', shell_quote($dst));
                die "$name: could not gzip $entry\n" if $st != 0;
            }
            else {
                copy($src, File::Spec->catfile($scratch, $entry)) or die "copy $src: $!\n";
            }
        }
        closedir $dh;
    }

    ### Staged under its own basename rather than a fixed name, because vcf_v7 has to be called
    ### mgp_REL2005_snps_indels.vcf.gz to reach the filename detection at :1435.
    my $gzip = -e File::Spec->catfile($dir, 'gzip_vcf');
    for my $vcf (sort glob File::Spec->catfile($dir, '*.vcf')) {
        my $base = basename($vcf);
        if ($gzip) {
            my $dst = File::Spec->catfile($scratch, "$base.gz");
            my $st  = system(join ' ', 'gzip', '-c', shell_quote($vcf), '>', shell_quote($dst));
            die "$name: could not gzip $base\n" if $st != 0;
        }
        else {
            copy($vcf, File::Spec->catfile($scratch, $base)) or die "copy $vcf: $!\n";
        }
    }
}

### readdir rather than glob, because empty_genome_folder ships only a .gitkeep - git cannot store an
### empty directory, and glob would skip the placeholder and stage nothing.
sub copy_tree {
    my ($src, $dst) = @_;
    make_path($dst);
    opendir my $dh, $src or die "opendir $src: $!\n";
    for my $entry (grep { $_ ne '.' && $_ ne '..' } readdir $dh) {
        my $from = File::Spec->catfile($src, $entry);
        my $to   = File::Spec->catfile($dst, $entry);
        if (-d $from) { copy_tree($from, $to) }
        else          { copy($from, $to) or die "copy $from: $!\n" }
    }
    closedir $dh;
}

###############################################################################
### Normalisation
###############################################################################

sub normalise {
    my ($name, $dir, $scratch, $exit, $impl_name) = @_;

    my @excludes = read_lines(File::Spec->catfile($dir, 'exclude'));
    my @present  = collect_files($scratch);

    my %out;
    ### Only whether the run aborted, not the number. Perl's die exits with $!, which reflects the last
    ### failed syscall: the same abort gave 2 here and 25 on another machine.
    $out{'exit_status'} = $exit == 0 ? "0\n" : "nonzero\n";
    ### manifest does more work here than in the SNPsplit suite: the four genome modes differ mainly in
    ### which of six output directories exist, and nothing else asserts that.
    $out{'manifest'} = join('', map { "$_\n" } @present);

    for my $rel (@present) {
        next if $rel eq 'stdout.log' || $rel eq 'stderr.log';
        next if excluded($rel, \@excludes);

        my $abs = File::Spec->catfile($scratch, $rel);

        ### Recorded as text under a .txt suffix, following the .bam -> .sam.txt precedent: the suffix
        ### signals that the file was transformed, and manifest still pins the archive's existence.
        my ($key, $text) = $rel =~ /\.gz$/
            ? ("$rel.txt", gunzip_to_text($abs))
            : ($rel, mask_paths(slurp($abs), $scratch, $impl_name));

        for my $rule (@NORMALISE) {
            next unless $key =~ $rule->{pattern};
            $text = mask_id_column($text) if $rule->{mask_id};
            $text = sort_lines($text)     if $rule->{sort_lines};
        }
        $out{$key} = $text;
    }

    my $stderr = mask_paths(slurp(File::Spec->catfile($scratch, 'stderr.log')), $scratch, $impl_name);
    $out{'run.log'} = sort_strain_lists($stderr);

    my $stdout = slurp(File::Spec->catfile($scratch, 'stdout.log'));
    if (length $stdout) {
        $stdout = mask_paths($stdout, $scratch, $impl_name);
        $stdout = sort_chrom_list($stdout);
        $stdout = sort_strain_lists($stdout);
        $out{'stdout.log'} = $stdout;
    }

    return %out;
}

### A zero-byte archive is the recorded result of poison_gzip, so an unreadable one becomes a
### comparable sentinel rather than aborting the suite.
sub gunzip_to_text {
    my ($abs) = @_;
    return "<unreadable>\n" if -z $abs;
    my $text = `gunzip -c @{[ shell_quote($abs) ]} 2>/dev/null`;
    return "<unreadable>\n" if $? != 0 || !length $text;
    return $text;
}

sub sort_lines {
    my ($text) = @_;
    return $text unless length $text;
    my @lines = split /^/, $text;
    ### A final line without a newline would sort into the middle of the file, so it is chomped with
    ### the rest and the newline handed back afterwards. The result is the input's lines, reordered.
    my $had_newline = $lines[-1] =~ /\n\z/ ? 1 : 0;
    chomp @lines;
    my $out = join "\n", sort @lines;
    $out .= "\n" if $had_newline;
    return $out;
}

### Field 0 is either a bare counter (:465, :499) or ${strain}_${counter} (:538). Only the digits vary
### between runs; the prefix is the only signal that a line came from the new-reference branch, and it
### is a documented output column, so it stays asserted.
sub mask_id_column {
    my ($text) = @_;
    my @out;
    for my $line (split /^/, $text) {
        $line =~ s/^([^\t]*?_)?(\d+)\t/(defined $1 ? $1 : '') . "<n>\t"/e;
        push @out, $line;
    }
    return join '', @out;
}

### detect_chroms returns keys %chrom unsorted (:1232), and the list is printed to stdout as one
### tab-separated line.
sub sort_chrom_list {
    my ($text) = @_;
    my @lines = split /^/, $text;
    for my $i (0 .. $#lines) {
        next unless $lines[$i] =~ /^Using the following chromosomes/;
        for my $j ($i + 1 .. $#lines) {
            next unless $lines[$j] =~ /\S/;
            chomp(my $list = $lines[$j]);
            $lines[$j] = join("\t", sort split /\t/, $list) . "\n";
            last;
        }
        last;
    }
    return join '', @lines;
}

### process_commandline prints the available strain names with `foreach my $strain (keys %strains)` at
### :1472, :1490, :1501, :1530 and :1541 - five paths, all unsorted, all to stdout. Keyed on the shape
### of the lines rather than on the fixture, so all five are covered.
###
### A run of one line sorts to itself, so a single-chromosome fixture's chromosome list - the only other
### bare-token stdout line - is unaffected.
sub sort_strain_lists {
    my ($text) = @_;
    my @lines = split /^/, $text;
    my @out;
    my @run;
    for my $line (@lines) {
        if ($line =~ /^[A-Za-z0-9_.-]+\n$/) { push @run, $line; next }
        push @out, sort @run if @run;
        @run = ();
        push @out, $line;
    }
    push @out, sort @run if @run;
    return join '', @out;
}

### Could not be copied verbatim: run_tests.pl's mask_paths closes over a file-scoped $samtools and
### carries a `Samtools path:` rule. This one takes the implementation name instead, to mask the line
### numbers in perl's own diagnostics.
sub mask_paths {
    my ($text, $scratch, $impl_name) = @_;
    my $real = defined $scratch ? abs_path($scratch) : undef;
    for my $p (grep { defined && length } ($real, $scratch)) {
        $text =~ s/\Q$p\E/<scratch>/g;
    }
    ### Four fixtures pin a die or warn that perl stamps with a line number, because those dies lack a
    ### trailing newline. Masking it keeps them from breaking on every edit above the reported line.
    ### The `<IN> line N` half is the input line number and is real behaviour, so it survives.
    $text =~ s/(\Q$impl_name\E line )\d+/$1<n>/g if defined $impl_name;
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

### Could not be copied verbatim: run_tests.pl's guard inspects alignment output. The genome analogue
### has to be scoped to the directories the tool generates, because the fixture's own staged genome/*.fa
### is also a .fa carrying sequence and would satisfy an unscoped check in every fixture, including the
### ones that produce nothing.
sub write_expected {
    my ($name, $dir, $actual, $exit) = @_;

    my $allow_empty = -e File::Spec->catfile($dir, 'allow_empty');

    if ($exit == 0 && !$allow_empty && generated_genome_is_empty($actual)) {
        warn "FAIL $name: the run succeeded but generated no genome sequence; refusing to record it.\n"
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
    printf "UPDATED %-24s %d expected files\n", $name, scalar keys %$actual;
    return 1;
}

sub generated_genome_is_empty {
    my ($actual) = @_;
    my @generated = grep { m{^[^/]*_(?:N-masked|full_sequence)/.*\.fa$} } keys %$actual;
    return 1 unless @generated;
    for my $key (@generated) {
        return 0 if grep { !/^>/ && /\S/ } split /^/, $actual->{$key};
    }
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

### Could not be copied verbatim: run_tests.pl's all_fixtures reads a file-scoped $fixture_dir.
sub all_fixtures {
    my ($fixtures) = @_;
    return sort map { basename($_) } grep { -d } glob File::Spec->catfile($fixtures, '*');
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

### expected/ is one flat directory, so nested output paths are encoded.
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
Differential regression harness for SNPsplit_genome_preparation.

  test/run_genome_tests.pl [OPTIONS] [FIXTURE ...]

  --update          Regenerate expected/ for the selected fixtures instead of diffing.
  --impl PATH       Implementation under test. Default: the repo's ./SNPsplit_genome_preparation.
  --keep            Leave scratch directories in place for inspection.
  -v, --verbose     Echo each command as it runs.
  -h, --help        This message.

  With no FIXTURE arguments, every directory under test/genome_fixtures/ runs.
  Exit status is 0 when all selected fixtures match, 1 otherwise.

  Requires gzip and gunzip on PATH. samtools is not used.

  Never run --update without reading the resulting diff: it will happily pin a bug.
USAGE
    exit $code;
}
