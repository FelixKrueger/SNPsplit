#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname basename);
use File::Spec ();
use Cwd qw(abs_path getcwd);

### Runs two versions of SNPsplit over the same real BAM and compares what they produced.
###
### The fixture suite proves that expected output did not change; this proves it on data the fixtures
### cannot imitate - millions of reads, every SNP shape, both strands. It is the check that caught
### nothing between 0.6.0 and 0.7.0 and is the only thing that could have contradicted #94.
###
###   test/bin/compare_versions.pl --bam sample.bam --snps all_SNPs_STRAIN_GRCm39.txt.gz \
###       --old 0.8.0 --new dev [--args '--paired --bisulfite']
###
### --old and --new are anything `git show` understands: a tag, a branch, a commit.

my ($bam,$snps,$old,$new,$extra,$outdir,$samtools,$keep,$help);
$old = '0.8.0';
$new = 'dev';
$extra = '';

GetOptions(
    'bam=s'      => \$bam,
    'snps=s'     => \$snps,
    'old=s'      => \$old,
    'new=s'      => \$new,
    'args=s'     => \$extra,
    'outdir=s'   => \$outdir,
    'samtools=s' => \$samtools,
    'keep'       => \$keep,
    'h|help'     => \$help,
) or usage(1);
usage(0) if $help;

usage(1) unless defined $bam and defined $snps;
$bam  = abs_path($bam);
$snps = abs_path($snps);
die "No such BAM: $bam\n"       unless -e $bam;
die "No such SNP file: $snps\n" unless -e $snps;

my $repo = abs_path(File::Spec->catdir(dirname(__FILE__), '..', '..'));

### SNPsplit tests --samtools_path with -e, so a bare program name is rejected. Resolved here rather
### than passed through, so the same binary is used for both runs and for the comparison.
$samtools = defined $samtools ? abs_path($samtools) : which('samtools');
die "samtools not found; pass --samtools PATH\n" unless defined $samtools && -x $samtools;
$outdir = defined $outdir ? abs_path($outdir)
                          : File::Spec->catdir(($ENV{TMPDIR} || '/tmp'), "snpsplit_compare_$$");

warn "Comparing '$old' against '$new'\n";
warn "BAM:       $bam\n";
warn "SNP file:  $snps\n";
warn "Workspace: $outdir\n\n";

my %ran;
for my $rev ($old, $new) {
    my $label = label_for($rev);
    my $dir   = File::Spec->catdir($outdir, $label);
    make_path($dir);

    ### All three scripts, because SNPsplit calls $RealBin/tag2sort with no override: mixing versions
    ### across that boundary would compare something nobody runs.
    for my $script (qw(SNPsplit tag2sort SNPsplit_genome_preparation)) {
        my $dst = File::Spec->catfile($dir, $script);
        run("git -C '$repo' show '$rev:$script' > '$dst'");
        chmod 0755, $dst;
    }

    my $work = File::Spec->catdir($dir, 'run');
    make_path($work);

    my $cwd = getcwd();
    chdir $work or die "chdir $work: $!\n";
    my $cmd = sprintf("'%s' --snp_file '%s' --samtools_path '%s' %s '%s' > stdout.log 2> stderr.log",
                      File::Spec->catfile($dir, 'SNPsplit'), $snps, $samtools, $extra, $bam);
    warn "[$label] $cmd\n";
    my $status = system($cmd);
    chdir $cwd or die "chdir back: $!\n";

    warn "[$label] exit status: @{[ $status >> 8 ]}\n\n";
    $ran{$label} = { dir => $work, exit => $status >> 8 };
}

my ($a, $b) = map { label_for($_) } ($old, $new);

print "\n", '=' x 72, "\n";
printf "exit status:  %s %d   |   %s %d\n", $a, $ran{$a}{exit}, $b, $ran{$b}{exit};
print '=' x 72, "\n\n";

### A run that produced nothing would otherwise sail through the comparison below with no files to
### compare and be reported as no change - the same vacuous pass this whole exercise exists to avoid.
for my $label ($a, $b) {
    next unless $ran{$label}{exit};
    print "$label exited $ran{$label}{exit}. Nothing was compared.\n";
    print "  see $ran{$label}{dir}/stderr.log\n";
    print "Workspace kept at $outdir\n";
    exit 2;
}

my $produced = keys %{ listing($ran{$a}{dir}) };
unless ($produced) {
    print "$a produced no output files. Nothing was compared.\n";
    print "Workspace kept at $outdir\n";
    exit 2;
}

### Alignment records are the thing that must not change. The header is compared separately because
### @PG records the argv, and 0.9.0 writes BAM with `samtools view -o` instead of a shell redirect, so
### the CL: field legitimately differs.
my $verdict = 0;
for my $rel (sort keys %{ { map { %{ listing($ran{$_}{dir}) } } ($a, $b) } }) {
    my $fa = File::Spec->catfile($ran{$a}{dir}, $rel);
    my $fb = File::Spec->catfile($ran{$b}{dir}, $rel);

    unless (-e $fa and -e $fb) {
        printf "%-46s %s\n", $rel, (-e $fa ? "only in $a" : "only in $b");
        $verdict = 1;
        next;
    }

    if ($rel =~ /\.bam$/) {
        my $records = same(qq('$samtools' view '$fa'), qq('$samtools' view '$fb'));
        my $header  = same(qq('$samtools' view -H '$fa' | grep -v '^\@PG'),
                           qq('$samtools' view -H '$fb' | grep -v '^\@PG'));
        printf "%-46s records %s   header (minus \@PG) %s\n", $rel,
               $records ? 'identical' : 'DIFFER', $header ? 'identical' : 'DIFFER';
        $verdict = 1 unless $records and $header;
    }
    elsif ($rel =~ /\.(txt|yaml|log)$/) {
        ### version, date_run and command differ by construction, and are the same three fields the
        ### fixture harness masks in the YAML.
        my $filter = q(grep -v -e 'ersion' -e 'date_run' -e 'command' -e ') . $outdir . q(');
        my $ok = same(qq($filter '$fa'), qq($filter '$fb'));
        printf "%-46s %s\n", $rel, $ok ? 'identical' : 'differs';
        $verdict = 1 if !$ok and $rel =~ /report\.txt$/;
    }
}

print "\n";
if ($verdict) {
    print "SOMETHING CHANGED. Read the differing files before releasing.\n";
}
else {
    print "No change in alignment records or reports. This is the expected result:\n";
    print "neither #99 nor #116 touches allele assignment.\n";
}
print "Workspace kept at $outdir\n" if $keep;
remove_tree($outdir) unless $keep;
exit $verdict;

###############################################################################

sub label_for { my $r = shift; $r =~ s/[^A-Za-z0-9._-]/_/g; return $r }

sub listing {
    my ($dir) = @_;
    my %f;
    for my $p (glob File::Spec->catfile($dir, '*')) {
        next unless -f $p;
        next if basename($p) =~ /^std(out|err)\.log$/;
        $f{ basename($p) } = 1;
    }
    return \%f;
}

sub same {
    my ($cmd_a, $cmd_b) = @_;
    my $out_a = `$cmd_a 2>/dev/null`;
    my $out_b = `$cmd_b 2>/dev/null`;
    return $out_a eq $out_b;
}

sub run {
    my ($cmd) = @_;
    system($cmd) == 0 or die "failed: $cmd\n";
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
Run two versions of SNPsplit over the same BAM and compare the output.

  test/bin/compare_versions.pl --bam FILE --snps FILE [OPTIONS]

  --bam FILE        Real alignment file to run both versions over. Required.
  --snps FILE       SNP annotation, as passed to --snp_file. Required.
  --old REV         Version to compare against. Default: 0.8.0
  --new REV         Version under test. Default: dev
  --args 'STRING'   Extra SNPsplit arguments, e.g. '--paired --bisulfite'
  --samtools PATH   samtools to use. Default: samtools on PATH
  --outdir DIR      Where to work. Default: a temporary directory
  --keep            Keep the workspace for inspection
  -h, --help        This message

  Exit status is 0 when nothing that matters changed.

  Alignment records are compared in full. BAM headers are compared with @PG lines
  removed, because @PG records the argv and 0.9.0 writes BAM through
  `samtools view -o` rather than a shell redirect.
USAGE
    exit $code;
}
