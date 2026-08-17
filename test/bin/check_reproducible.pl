#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname basename);
use File::Spec ();
use File::Find ();
use Cwd qw(abs_path getcwd);

### Asserts that SNPsplit_genome_preparation produces byte-identical output on identical input.
###
### The fixture suite cannot do this. Six outputs were written straight out of Perl hashes (#104), so
### the runner sorts them before comparison - which means the suite passes whether or not the source
### sorts them, and a regression that removed the sorting would surface as nothing at all here and as
### an unreproducible genome build for the user.
###
### So this runs the tool twice under two fixed hash seeds and compares the raw output. Sorting is kept
### in the runner as well, because --impl has to stay usable against an implementation that has its own
### ideas about iteration order.
###
###   test/bin/check_reproducible.pl [--impl PATH]

my $test_root = abs_path(File::Spec->catdir(dirname(__FILE__), '..'));
my $repo_root = abs_path(File::Spec->catdir($test_root, '..'));

my $impl = File::Spec->catfile($repo_root, 'SNPsplit_genome_preparation');
if (@ARGV >= 2 && $ARGV[0] eq '--impl') { $impl = abs_path($ARGV[1]) }
die "implementation not found or not executable: $impl\n" unless -x $impl;

### dual_hybrid produces every output that was affected, including the annotation file whose ID column
### was nondeterministic in content rather than just in order.
my $fixture = File::Spec->catdir($test_root, 'genome_fixtures', 'dual_hybrid');
die "fixture not found: $fixture\n" unless -d $fixture;

my $work = File::Spec->catdir(($ENV{TMPDIR} || '/tmp'), "gp_reproducible_$$");
remove_tree($work) if -d $work;

my @seeds = (1, 4242);
my %tree;

for my $seed (@seeds) {
    my $dir = File::Spec->catdir($work, "seed_$seed");
    make_path(File::Spec->catdir($dir, 'genome'));
    copy_in(File::Spec->catdir($fixture, 'genome'), File::Spec->catdir($dir, 'genome'));
    copy_in($fixture, $dir, qr/\.vcf$/);

    my $args = read_args(File::Spec->catfile($fixture, 'args'));
    my $cwd  = getcwd();
    chdir $dir or die "chdir $dir: $!\n";
    local $ENV{PERL_HASH_SEED}    = $seed;
    local $ENV{PERL_PERTURB_KEYS} = 2;
    local $ENV{SNPSPLIT_NO_SLEEP} = 1;
    my $status = system("'$impl' $args > stdout.log 2> stderr.log");
    chdir $cwd or die "chdir back: $!\n";
    die "the run failed under seed $seed (exit @{[ $status >> 8 ]}); see $dir/stderr.log\n"
        if $status != 0;

    $tree{$seed} = collect($dir);
}

### The logs carry the scratch path and so always differ; the staged input is not output.
my @differing;
my @keys = sort keys %{ $tree{ $seeds[0] } };
for my $rel (@keys) {
    next if $rel =~ m{^(?:stdout\.log|stderr\.log|genome/)} || $rel =~ /\.vcf$/;
    my $a = $tree{ $seeds[0] }{$rel};
    my $b = $tree{ $seeds[1] }{$rel};
    unless (defined $b) { push @differing, "$rel (missing under seed $seeds[1])"; next }
    push @differing, $rel unless $a eq $b;
}

printf "Compared %d output files across hash seeds %s\n", scalar @keys, join(' and ', @seeds);

if (@differing) {
    print "NOT REPRODUCIBLE - these differ between two runs on identical input:\n";
    print "  $_\n" for @differing;
    print "\nA hash-iteration site is unsorted. See #104.\n";
    print "Kept for inspection: $work\n";
    exit 1;
}

remove_tree($work);
print "Reproducible: every output file is byte-identical under both seeds.\n";
exit 0;

###############################################################################

sub collect {
    my ($dir) = @_;
    my %out;
    File::Find::find(
        {
            no_chdir => 1,
            wanted   => sub {
                my $path = $File::Find::name;
                return unless -f $path;
                my $rel = File::Spec->abs2rel($path, $dir);
                ### Compared as text, so a gzip header's timestamp cannot mask the content.
                $out{$rel} = $rel =~ /\.gz$/ ? `gunzip -c '$path'` : slurp($path);
            },
        },
        $dir
    );
    return \%out;
}

sub copy_in {
    my ($from, $to, $only) = @_;
    opendir my $dh, $from or die "opendir $from: $!\n";
    for my $entry (grep { $_ ne '.' && $_ ne '..' } readdir $dh) {
        my $src = File::Spec->catfile($from, $entry);
        next unless -f $src;
        next if $only && $entry !~ $only;
        open my $in,  '<', $src or die "read $src: $!\n";
        open my $out, '>', File::Spec->catfile($to, $entry) or die "write $entry: $!\n";
        print {$out} do { local $/; <$in> };
        close $in;
        close $out;
    }
    closedir $dh;
}

sub read_args {
    my ($file) = @_;
    open my $fh, '<', $file or die "read $file: $!\n";
    my @lines = grep { /\S/ && !/^\s*#/ } map { chomp; $_ } <$fh>;
    close $fh;
    return join ' ', @lines;
}

sub slurp {
    my ($file) = @_;
    open my $fh, '<', $file or return '';
    local $/;
    my $c = <$fh>;
    close $fh;
    return defined $c ? $c : '';
}
