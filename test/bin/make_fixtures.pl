#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use File::Spec ();
use File::Basename qw(dirname);
use Cwd qw(abs_path getcwd);

### Authors every fixture: reference, SNP annotation, SAM records and markers.
###
### MD tags are derived with `samtools calmd` against the committed N-masked reference rather than
### computed by hand, because nothing downstream validates an MD tag - a self-consistent tag pointing
### at the wrong coordinate produces plausible-looking expected output that pins wrong behaviour.
###
### Re-run after editing the read tables, then re-record the expected output and read every diff:
###   test/bin/make_fixtures.pl && test/run_tests.pl --update

my $test_root = abs_path(File::Spec->catdir(dirname(__FILE__), '..'));
my $fixtures  = File::Spec->catdir($test_root, 'fixtures');
my $samtools  = $ENV{SAMTOOLS} || 'samtools';

### --- Shared reference ------------------------------------------------------------------------
###
### SNP positions are masked to N exactly as SNPsplit_genome_preparation would, because the whole
### mechanism depends on the aligner reporting N as a mismatch in the MD tag.
###
### Spacing is deliberate - a 50 bp read starting at each of these positions spans exactly the SNPs
### named, and nothing else:
###   pos   6 (6-55)    -> SNP 30 only          pos 110 (110-159) -> SNP 130 only
###   pos  20 (20-69)   -> SNP 30 and SNP 60    pos 180 (180-229) -> SNP 200 only
###   pos  40 (40-89)   -> SNP 60 only          pos 230 (230-279) -> no SNP at all
###                                             pos 260 (260-309) -> the unlisted N at 280

my %STD = (
    len      => 320,
    snps     => { 30 => ['A','G'], 60 => ['T','C'], 130 => ['G','T'], 200 => ['C','A'] },
    ### Masked in the reference but absent from snps.txt: the only way to reach the
    ### "N was not present in the list of known SNPs" counter.
    unlisted => { 280 => 'G' },
);

### Bisulfite needs one SNP of every shape score_bisulfite_SNPs dispatches on, spaced 60 apart so a
### 50 bp read isolates each:
###   30 C/T  -> ref eq 'C'                 330 A/T  -> touches neither C nor G
###   90 C/G  -> ref eq 'C', snp eq 'G'     390 G/C  -> snp eq 'C', ref eq 'G'
###  150 T/C  -> snp eq 'C'                 450 G/T  -> ref eq 'G', bottom-strand arm
###  210 G/A  -> ref eq 'G'                 510 T/G  -> snp eq 'G', bottom-strand arm
###  270 A/G  -> snp eq 'G'
my %BS = (
    len  => 540,
    snps => {
        30  => ['C','T'], 90  => ['C','G'], 150 => ['T','C'],
        210 => ['G','A'], 270 => ['A','G'], 330 => ['A','T'],
        390 => ['G','C'], 450 => ['G','T'], 510 => ['T','G'],
    },
    unlisted => {},
);

my $PG_BOWTIE  = qq(\@PG\tID:bowtie2\tPN:bowtie2\tVN:2.4.5\tCL:"synthetic fixture");
### Needs ID:Bismark for $isBismark, plus whitespace-delimited -1 and -2 for pairing inference.
my $PG_BISMARK = qq(\@PG\tID:Bismark\tPN:Bismark\tVN:0.24.0\tCL:"bismark --genome /g -1 r1.fq -2 r2.fq");

my $OT   = 'XR:Z:CT	XG:Z:CT';
my $OB   = 'XR:Z:CT	XG:Z:GA';
my $CTOT = 'XR:Z:GA	XG:Z:CT';
my $CTOB = 'XR:Z:GA	XG:Z:GA';

### --- Fixtures -------------------------------------------------------------------------------

my %F;

##### Tier 1 ###################################################################################

$F{se_basic} = {
    args  => '--single_end',
    notes => 'ref/alt/neither base at a SNP, multi-SNP agreement and conflict, no-SNP read, unlisted N',
    reads => [
        { q => 'r_g1',         pos => 6,   subs => { 30 => 'A' } },                # ref base      -> G1
        { q => 'r_g2',         pos => 6,   subs => { 30 => 'G' } },                # alt base      -> G2
        { q => 'r_neither',    pos => 6,   subs => { 30 => 'C' } },                # third base    -> UA, no_snp
        { q => 'r_no_snp',     pos => 230 },                                        # no N in MD    -> UA
        { q => 'r_two_g1',     pos => 20,  subs => { 30 => 'A', 60 => 'T' } },     # both ref      -> G1
        { q => 'r_conflict',   pos => 20,  subs => { 30 => 'A', 60 => 'C' } },     # ref + alt     -> CF
        { q => 'r_unlisted_n', pos => 260 },                                        # N not in list -> no_snp_found
        { q => 'r_snp3_g2',    pos => 110, subs => { 130 => 'T' } },               # alt at SNP3   -> G2
        { q => 'r_snp4_g1',    pos => 180, subs => { 200 => 'C' } },               # ref at SNP4   -> G1
    ],
};

$F{se_skipped} = {
    args     => '--single_end',
    notes    => 'unmapped, hard-clipped and missing-MD reads are all skipped before allele assignment',
    strip_md => ['r_nomd'],
    reads    => [
        { q => 'r_ok',       pos => 6, subs => { 30 => 'A' } },                     # keeps output non-empty
        { q => 'r_unmapped', pos => 6, flag => 4, subs => { 30 => 'A' } },          # 0x4 -> unmapped counter
        { q => 'r_hardclip', pos => 6, cigar => '5H45M', subs => { 30 => 'A' } },   # H -> hardclipped (stderr only)
        { q => 'r_nomd',     pos => 6, subs => { 30 => 'A' } },                     # MD stripped -> counted nowhere
    ],
};

### princess_paired_end collapses G1/G1, G1/UA and UA/G1 into one counter, so unlike Hi-C the
### ordering variants are coverage-equivalent here and one of each suffices.
### Singletons sit in the middle so a differently-named read follows each one and the routing block
### inside the loop runs; a singleton placed last is handled by the end-of-file block instead, so one
### is kept there too. Both conflicting arms are covered: a G1/G2 pair and a pair with a CF mate.
$F{pe_basic} = {
    args  => '--paired --no_sorting --conflicting',
    notes => 'all four pair outcomes, both conflicting arms, and all four singleton routes in both blocks',
    reads => [
        pe_pair('p_g1_g1', 'G1', 'G1'),
        { q => 's_g1', flag => 65, %{ read_for('G1') } },
        { q => 's_g2', flag => 65, %{ read_for('G2') } },
        { q => 's_ua', flag => 65, %{ read_for('UA') } },
        { q => 's_cf', flag => 65, %{ read_for('CF') } },
        pe_pair('p_g1_g2', 'G1', 'G2'),
        pe_pair('p_cf_ua', 'CF', 'UA'),
        pe_pair('p_ua_ua', 'UA', 'UA'),
        pe_pair('p_g1_ua', 'G1', 'UA'),
        { q => 's_last', flag => 65, %{ read_for('G1') } },
    ],
};

### All ten branches of process_HiC_paired_end. It counts G1/UA separately from UA/G1 (likewise
### G2/UA vs UA/G2 and G1/G2 vs G2/G1) behind shared output files, so both orderings are required or
### three branches stay unreached while the expected output still looks complete. CF is tested before the
### mixed cases, so CF/UA reaches the conflicting branch.
$F{hic} = {
    args  => '--hic --conflicting',
    notes => 'all ten pair branches of process_HiC_paired_end, both orderings of every mixed pair',
    reads => [ map { hic_pair(@$_) }
        ['h_ua_ua','UA','UA'], ['h_g1_g1','G1','G1'], ['h_g2_g2','G2','G2'],
        ['h_cf_ua','CF','UA'], ['h_g1_ua','G1','UA'], ['h_ua_g1','UA','G1'],
        ['h_g2_ua','G2','UA'], ['h_ua_g2','UA','G2'], ['h_g1_g2','G1','G2'],
        ['h_g2_g1','G2','G1'],
    ],
};

$F{samtools_path} = {
    args   => '--single_end --samtools_path <REAL_SAMTOOLS>',
    notes  => 'regression test: every samtools call must honour --samtools_path, reads included',
    poison => 1,
    reads  => $F{se_basic}{reads},
};

##### Tier 2 ###################################################################################

### The %insertions_accounted_for / %softclips_accounted_for arms are reachable only with two or more
### Ns and an insertion or soft-clip before both, because STEP II re-walks the CIGAR for every N.
### The last two reads put the masked position *after* the gap. With the N inside the first M block the
### walker stops before it ever reaches the D or N operation, so the adjustment that maps a
### post-junction SNP onto its genomic coordinate never runs.
$F{cigar_complex} = {
    args  => '--single_end',
    notes => 'complex CIGARs across all three assignment outcomes, plus masked positions after a splice and a deletion',
    reads => [
        { q => 'r_ins',        pos => 6,   cigar => '25M3I22M',   subs => { 30 => 'A' } },   # ref   -> G1
        { q => 'r_ins_alt',    pos => 6,   cigar => '25M3I22M',   subs => { 30 => 'G' } },   # alt   -> G2
        { q => 'r_ins_third',  pos => 6,   cigar => '25M3I22M',   subs => { 30 => 'C' } },   # third -> unassignable
        { q => 'r_ins_unl',    pos => 256, cigar => '25M3I22M' },                            # unlisted N at 280
        { q => 'r_del',        pos => 6,   cigar => '30M2D20M',   subs => { 30 => 'A' } },
        { q => 'r_soft_lead',  pos => 6,   cigar => '5S45M',      subs => { 30 => 'A' } },
        { q => 'r_soft_tail',  pos => 6,   cigar => '45M5S',      subs => { 30 => 'A' } },
        { q => 'r_splice',     pos => 6,   cigar => '25M100N25M', subs => { 30 => 'A' } },
        { q => 'r_two_n_si',   pos => 20,  cigar => '4S6M3I35M',  subs => { 30 => 'A', 60 => 'T' } },
        { q => 'r_after_gap',  pos => 65,  cigar => '25M20N25M',  subs => { 130 => 'T' } },  # N in the second exon
        { q => 'r_after_del',  pos => 100, cigar => '20M5D30M',   subs => { 130 => 'T' } },  # N after the deletion
    ],
};

### A D operation over an N position puts ^N in the MD tag; SNPsplit counts those deletions and calls
### the read conflicting.
$F{n_deletions} = {
    args  => '--single_end',
    notes => 'one deleted N position, then two in the same read (multi_N_deletion)',
    reads => [
        { q => 'r_ok',        pos => 6,  subs => { 30 => 'A' } },
        { q => 'r_del_n_1',   pos => 6,  cigar => '24M1D25M' },          # MD 24^N25    -> CF
        { q => 'r_del_n_2',   pos => 20, cigar => '10M1D29M1D9M' },      # MD 10^N29^N9 -> multi
    ],
};

### A C>T SNP is unusable on the top strands, where conversion makes C and T indistinguishable, but
### fully usable on the bottom strands, which read the complement. That asymmetry is why rejection is
### per-read rather than at load time.
###
### Reads carry the reference base, the SNP base and a third base on both strands, because a fixture
### where every read matches the reference reaches only the genome-1 leaves and leaves the whole
### genome-2 half of the subroutine untested.
$F{bisulfite_se} = {
    args  => '--bisulfite --single_end',
    notes => 'every SNP shape and assignment leaf of score_bisulfite_SNPs, across all four Bismark strands',
    ref   => \%BS,
    reads => [ bs_reads() ],
};

$F{pe_singletons} = {
    args  => '--paired --no_sorting --conflicting --singletons',
    notes => 'same pairs as pe_basic, routed into the separate _st singleton output files',
    reads => $F{pe_basic}{reads},
};

### Name-sorting is gated on --paired without --no_sorting; it needs neither Bismark nor bisulfite.
### The samtools-written intermediate is excluded from content comparison because its ordering and
### @HD rewriting are samtools' behaviour, not SNPsplit's.
$F{pe_namesort} = {
    args    => '--paired --conflicting',
    notes   => 'mates deliberately non-adjacent, so sort_by_name_paired_end runs and renames outputs',
    exclude => ['*.sortedByName.bam'],
    reads   => [
        { q => 'n_a', flag => 65,  pos => 6,   subs => { 30 => 'A' } },
        { q => 'n_b', flag => 65,  pos => 6,   subs => { 30 => 'G' } },
        { q => 'n_a', flag => 129, pos => 6,   subs => { 30 => 'A' } },
        { q => 'n_b', flag => 129, pos => 6,   subs => { 30 => 'G' } },
    ],
};

### sam2bam_convert is triggered by the input filename, independent of --sam (which is broken end to
### end and deliberately not exercised here).
$F{sam_input} = {
    args     => '--single_end',
    notes    => 'SAM input converted by sam2bam_convert; --sam is a separate and broken mechanism',
    feed_sam => 1,
    exclude  => ['sam_input.bam'],
    reads    => $F{se_basic}{reads},
};

### --sam has to work end to end: SAM in, SAM out, sorting report merged, no intermediate left behind.
$F{sam_output} = {
    args  => '--single_end --sam',
    notes => 'SAM output throughout, including the sorting report merge and cleanup',
    reads => $F{se_basic}{reads},
};

### A .sam input with no --samtools_path, so the failing shim is the samtools that sam2bam_convert
### finds. Tests that a failed SAM to BAM conversion aborts rather than being reported as success.
$F{sam2bam_fails} = {
    args        => '--single_end',
    notes       => 'a failed SAM to BAM conversion aborts the run',
    feed_sam    => 1,
    poison      => 1,
    allow_empty => 1,
    reads       => $F{se_basic}{reads},
};

### The sorting stage is replaced with one that always fails, so this tests the abort itself rather
### than any particular reason for failing.
$F{tag2sort_fails} = {
    args         => '--single_end',
    notes        => 'a failing sorting stage aborts the run instead of being reported as success',
    break_sorter => 1,
    allow_empty  => 1,
    reads        => $F{se_basic}{reads},
};

### Auto-detection has two arms. A Bismark @PG without -1/-2 infers single-end and still sets
### --bisulfite, which no other fixture reaches.
$F{bismark_autodetect_se} = {
    args  => '',
    notes => 'ID:Bismark @PG without -1/-2 infers single-end and auto-sets --bisulfite',
    pg    => qq(\@PG\tID:Bismark\tPN:Bismark\tVN:0.24.0\tCL:"bismark --genome /g reads.fq"),
    reads => [ { q => 'a_single', pos => 6, subs => { 30 => 'A' }, tags => $OT } ],
};

$F{skip_tag2sort} = {
    args  => '--single_end --skip_tag2sort',
    notes => 'tagging only: SNPsplit exits before the sorting stage, so no genome*.bam appear',
    reads => $F{se_basic}{reads},
};

### H is filtered before the CIGAR walker, but = and X are not, so an --eqx aligner's output aborts
### the run. Pinned by exit status rather than by output.
$F{cigar_eqx} = {
    args        => '--single_end',
    notes       => '=/X CIGAR operations abort the CIGAR walker; pinned via exit_status',
    allow_empty => 1,
    reads       => [ { q => 'r_eqx', pos => 6, cigar => '24=1X25=', subs => { 30 => 'A' } } ],
};

$F{header_only} = {
    args        => '--single_end',
    notes       => 'zero alignments: the $count == 0 branch and its N/A percentage formatting',
    allow_empty => 1,
    reads       => [],
};

### Auto-detection sets --bisulfite and infers pairing from the @PG line, then test_positional_sorting
### sets --no_sort because the mates are adjacent. No flags are passed at all.
$F{bismark_autodetect} = {
    args  => '',
    notes => 'ID:Bismark @PG auto-sets --bisulfite, infers --paired from -1/-2, and sets --no_sort',
    pg    => $PG_BISMARK,
    reads => [
        { q => 'a_pair', flag => 65,  pos => 6, subs => { 30 => 'A' }, tags => $OT },
        { q => 'a_pair', flag => 129, pos => 6, subs => { 30 => 'A' }, tags => $OT },
    ],
};

### Two input files are the only way to reach "Using already stored SNP information", and they expose
### the counters that are never reset between files.
$F{multi_input} = {
    args   => '--single_end',
    notes  => 'second input file reuses the stored SNPs; pins the cross-file counter carry-over',
    reads  => [ { q => 'm1_g1', pos => 6, subs => { 30 => 'A' } },
                { q => 'm1_g2', pos => 6, subs => { 30 => 'G' } } ],
    reads2 => [ { q => 'm2_g1', pos => 6, subs => { 30 => 'A' } },
                { q => 'm2_ua', pos => 230 } ],
};

### With --output_dir the -e guards on the sorting report look in the wrong directory, so the merge is
### skipped and SNPsplit_sort.yaml survives. Pinned as current behaviour.
$F{output_dir} = {
    args  => '--single_end --output_dir out',
    notes => 'output_dir: sorting report is not merged and SNPsplit_sort.yaml is left behind',
    reads => $F{se_basic}{reads},
};

### With no G1/G1 pair, $genome1 is never initialised, so the report prints a blank count and the YAML
### an empty value. Pinned as current behaviour.
$F{hic_no_g1} = {
    args  => '--hic --conflicting',
    notes => 'Hi-C without a G1/G1 pair: blank genome 1 count in the report and an empty YAML value',
    reads => [ map { hic_pair(@$_) } ['n_ua_ua','UA','UA'], ['n_g2_g2','G2','G2'] ],
};

### --- Emit -----------------------------------------------------------------------------------

for my $name (sort keys %F) {
    my $spec = $F{$name};
    my $cfg  = $spec->{ref} || \%STD;
    my ($plain, $masked) = make_ref($cfg);

    my $dir = File::Spec->catdir($fixtures, $name);
    make_path($dir);

    my $fa = File::Spec->catfile($dir, 'chr1.fa');
    write_file($fa, fasta($masked));
    write_file(File::Spec->catfile($dir, 'snps.txt'), snp_table($cfg));
    write_file(File::Spec->catfile($dir, 'args'),   ($spec->{args} || '') . "\n");
    write_file(File::Spec->catfile($dir, 'README'), "$spec->{notes}\n");

    my %marker_key = (poison_shim => 'poison', break_tag2sort => 'break_sorter');
    for my $marker (qw(poison_shim feed_sam allow_empty break_tag2sort)) {
        my $f = File::Spec->catfile($dir, $marker);
        $spec->{ $marker_key{$marker} || $marker } ? write_file($f, '') : unlink $f;
    }
    my $exc = File::Spec->catfile($dir, 'exclude');
    $spec->{exclude} ? write_file($exc, join("\n", @{ $spec->{exclude} }) . "\n") : unlink $exc;

    system("$samtools faidx '$fa'") == 0 or die "faidx failed for $name\n";

    my $n = 0;
    for my $key (qw(reads reads2)) {
        next unless $spec->{$key};
        my $out = $key eq 'reads' ? 'input.sam' : 'input2.sam';
        my $raw = File::Spec->catfile($dir, ".raw.sam");
        write_file($raw, sam_text($spec->{$key}, $plain, $cfg, $spec->{pg} || $PG_BOWTIE));

        ### calmd records its own argv in the @PG CL: field, so it runs from inside the fixture
        ### directory with relative names: absolute paths there would put whoever generated the
        ### fixtures into all 19 committed input files.
        my $md = do {
            my $cwd = getcwd();
            chdir $dir or die "chdir $dir: $!\n";
            my $text = `$samtools calmd .raw.sam chr1.fa 2>/dev/null`;
            chdir $cwd or die "chdir back: $!\n";
            $text;
        };
        die "calmd failed for $name/$out\n" if $? != 0 || !length $md;
        $md = strip_md_tag($md, $spec->{strip_md}) if $spec->{strip_md};

        write_file(File::Spec->catfile($dir, $out), $md);
        unlink $raw;
        $n += scalar @{ $spec->{$key} };
    }
    unlink "$fa.fai";

    printf "%-20s %2d records  args: %s\n", $name, $n, ($spec->{args} || '(auto-detect)');
}

print "\nNow re-record the expected output and read every diff:  test/run_tests.pl --update\n";

###############################################################################

### Deterministic sequence; SNP positions are forced to their reference base so the annotation and
### the sequence cannot disagree.
sub make_ref {
    my ($cfg) = @_;
    my @bases = qw(A C G T);
    my $seq   = '';
    my $x     = 12345;
    for (1 .. $cfg->{len}) {
        $x = ($x * 1103515245 + 12345) % 2147483648;
        $seq .= $bases[ ($x >> 16) % 4 ];
    }
    substr($seq, $_ - 1, 1) = $cfg->{snps}{$_}[0]  for keys %{ $cfg->{snps} };
    substr($seq, $_ - 1, 1) = $cfg->{unlisted}{$_} for keys %{ $cfg->{unlisted} };

    my $masked = $seq;
    substr($masked, $_ - 1, 1) = 'N' for (keys %{ $cfg->{snps} }, keys %{ $cfg->{unlisted} });
    return ($seq, $masked);
}

sub fasta {
    my ($seq) = @_;
    my $out = ">chr1\n";
    $out .= substr($seq, $_, 60) . "\n" for map { $_ * 60 } 0 .. int((length($seq) - 1) / 60);
    return $out;
}

sub snp_table {
    my ($cfg) = @_;
    my ($out, $i) = ('', 0);
    for my $pos (sort { $a <=> $b } keys %{ $cfg->{snps} }) {
        $i++;
        my ($ref, $alt) = @{ $cfg->{snps}{$pos} };
        $out .= join("\t", $i, 'chr1', $pos, 1, "$ref/$alt") . "\n";
    }
    return $out;
}

### Every SNP shape gets a reference-base, SNP-base and third-base read on both the top and the bottom
### strand. Each read starts 24 bp before its SNP, which puts the masked position at read offset 24 and
### keeps the neighbouring SNPs (60 bp away) outside the 50 bp span.
sub bs_reads {
    my @reads;
    for my $p (sort { $a <=> $b } keys %{ $BS{snps} }) {
        my ($ref, $alt) = @{ $BS{snps}{$p} };
        my $third = third_base($ref, $alt);
        my $start = $p - 24;
        for my $case (['ref', $ref], ['alt', $alt], ['third', $third]) {
            my ($label, $base) = @$case;
            push @reads,
                { q => "s${p}_ot_$label", pos => $start, subs => { $p => $base }, tags => $OT },
                { q => "s${p}_ob_$label", pos => $start, subs => { $p => $base }, tags => $OB };
        }
    }
    ### The remaining two Bismark strand values, and one indel read so the complex-CIGAR walker also
    ### dispatches into the bisulfite scorer.
    push @reads,
        { q => 'ctot_ct',  pos => 6, subs => { 30 => 'C' }, tags => $CTOT },
        { q => 'ctob_ct',  pos => 6, subs => { 30 => 'C' }, tags => $CTOB },
        { q => 'indel_ob', pos => 6, cigar => '25M3I22M', subs => { 30 => 'C' }, tags => $OB };
    return @reads;
}

sub third_base {
    my ($ref, $alt) = @_;
    for my $b (qw(A C G T)) { return $b if $b ne $ref && $b ne $alt }
}

sub pe_pair {
    my ($q, $t1, $t2) = @_;
    return ({ q => $q, flag => 65, %{ read_for($t1) } }, { q => $q, flag => 129, %{ read_for($t2) } });
}

sub hic_pair {
    my ($q, $t1, $t2) = @_;
    return ({ q => $q, flag => 65, %{ read_for($t1) } }, { q => $q, flag => 129, %{ read_for($t2) } });
}

sub read_for {
    my ($tag) = @_;
    return { pos => 6,   subs => { 30 => 'A' } }              if $tag eq 'G1';
    return { pos => 6,   subs => { 30 => 'G' } }              if $tag eq 'G2';
    return { pos => 230 }                                     if $tag eq 'UA';
    return { pos => 20,  subs => { 30 => 'A', 60 => 'C' } }   if $tag eq 'CF';
    die "unknown tag $tag\n";
}

sub sam_text {
    my ($reads, $plain, $cfg, $pg) = @_;

    ### The @PG line is required, not cosmetic: check_for_bs only reads --single_end inside its @PG
    ### branch, so a header without one dies before any allele tagging happens.
    my $out = "\@HD\tVN:1.6\tSO:unsorted\n\@SQ\tSN:chr1\tLN:$cfg->{len}\n$pg\n";

    for my $r (@$reads) {
        my $cigar = $r->{cigar} || '50M';
        my $seq   = build_seq($r->{pos}, $cigar, $r->{subs} || {}, $plain, $r->{q});
        my @f = ($r->{q}, ($r->{flag} || 0), 'chr1', $r->{pos}, 42, $cigar,
                 '*', 0, 0, $seq, 'I' x length($seq));
        push @f, $r->{tags} if $r->{tags};
        $out .= join("\t", @f) . "\n";
    }
    return $out;
}

### Walks the CIGAR so insertions and soft-clips contribute query bases the reference does not have,
### and deletions and splices advance the reference without contributing any.
sub build_seq {
    my ($pos, $cigar, $subs, $plain, $qname) = @_;
    my ($seq, $ref) = ('', $pos);
    my %aligned;
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        ### index() rather than a nested match: a regex inside the loop body would reset $1 and $2.
        my ($n, $op) = ($1, $2);
        if (index('M=X', $op) >= 0) {
            for my $i (0 .. $n - 1) {
                my $p = $ref + $i;
                $aligned{$p} = 1;
                $seq .= exists $subs->{$p} ? $subs->{$p} : substr($plain, $p - 1, 1);
            }
            $ref += $n;
        }
        elsif ($op eq 'I' or $op eq 'S') { $seq .= 'A' x $n }
        elsif ($op eq 'D' or $op eq 'N') { $ref += $n }
    }
    ### Every substituted position must be covered by an M/=/X operation, not merely inside the read's
    ### reference span: one that lands in a D or N gap contributes no base and would be lost in silence.
    for my $p (keys %$subs) {
        die "read $qname CIGAR $cigar does not align position $p\n" unless $aligned{$p};
    }
    return $seq;
}

sub strip_md_tag {
    my ($text, $qnames) = @_;
    my %want = map { $_ => 1 } @$qnames;
    my @out;
    for my $line (split /^/, $text) {
        if ($line !~ /^\@/) {
            my ($q) = split /\t/, $line;
            $line =~ s/\tMD:Z:\S+// if $want{$q};
        }
        push @out, $line;
    }
    return join '', @out;
}

sub write_file {
    my ($path, $content) = @_;
    open my $fh, '>', $path or die "write $path: $!\n";
    print {$fh} $content;
    close $fh;
}
