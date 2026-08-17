#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path remove_tree);
use File::Spec ();
use File::Basename qw(dirname);
use Cwd qw(abs_path);

### Authors every genome-preparation fixture: references, VCFs, pre-made SNP files, args and markers.
###
### Every VCF REF base is read out of the reference the VCF is written against, never typed, and
### asserted afterwards. A REF that disagrees with the reference is skipped silently by
### create_modified_chromosome:626 - the counter that records it is never reported - so a mistyped REF
### produces a fixture that records a plausible zero-SNP result and looks correct.
###
### Re-run after editing the tables, then re-record and read every diff:
###   test/bin/make_genome_fixtures.pl && test/run_genome_tests.pl --update

my $test_root = abs_path(File::Spec->catdir(dirname(__FILE__), '..'));
my $fixtures  = File::Spec->catdir($test_root, 'genome_fixtures');

### --- Reference sequence -----------------------------------------------------------------------
###
### A small LCG rather than a repeating motif, so a masked position is visible against its neighbours.
### The multiplier is kept small on purpose: the products stay well inside 2^31, so the sequence is
### identical on a 32-bit perl where a larger multiplier would drift into floating point.
sub sequence {
    my ($len, $seed) = @_;
    my @bases = qw(A C G T);
    my $x     = $seed;
    my $s     = '';
    for (1 .. $len) {
        $x = ($x * 75 + 74) % 65537;
        $s .= $bases[$x % 4];
    }
    return $s;
}

sub fasta {
    my ($name, $seq) = @_;
    my $out = ">$name\n";
    $out .= substr($seq, $_ * 60, 60) . "\n" for 0 .. int((length($seq) - 1) / 60);
    return $out;
}

### The ALT must differ from the REF or create_modified_chromosome:619 counts the position as already
### carrying the SNP and skips it.
sub alt_for {
    my ($ref, $offset) = @_;
    my @bases = grep { $_ ne $ref } qw(A C G T);
    return $bases[ ($offset || 0) % scalar @bases ];
}

my $SEQ  = sequence(320, 7);      # chromosome '1' in most fixtures
my $SEQ2 = sequence(280, 23);     # a second chromosome where one is needed
my $SEQ3 = sequence(200, 91);     # a third, deliberately absent from every VCF

### --- Fixture helpers --------------------------------------------------------------------------

my @written;

sub fixture {
    my (%f) = @_;
    my $dir = File::Spec->catdir($fixtures, $f{name});
    make_path($dir);

    ### Only the inputs this script owns are cleared. expected/ is deliberately left alone: wiping the
    ### fixture directory wholesale would delete every recorded output on any re-run, and the recording
    ### that followed would look like a clean pass while 28 fixtures had silently lost their baseline.
    remove_tree(File::Spec->catdir($dir, $_)) for qw(genome snps_in);
    unlink glob File::Spec->catfile($dir, '*.vcf');
    unlink grep { -e } map { File::Spec->catfile($dir, $_) }
        qw(args README pre_existing gzip_vcf poison_gzip allow_empty exclude);

    write_file(File::Spec->catfile($dir, 'args'),   "$f{args}\n");
    write_file(File::Spec->catfile($dir, 'README'), "$f{readme}\n");

    for my $marker (@{ $f{markers} || [] }) {
        write_file(File::Spec->catfile($dir, $marker), '');
    }
    if ($f{pre_existing}) {
        write_file(File::Spec->catfile($dir, 'pre_existing'), join('', map { "$_\n" } @{ $f{pre_existing} }));
    }

    ### genome/ is a hash of filename => contents; an empty hash still creates the directory, which is
    ### what empty_genome_folder needs.
    if ($f{genome}) {
        my $gdir = File::Spec->catdir($dir, 'genome');
        make_path($gdir);
        write_file(File::Spec->catfile($gdir, $_), $f{genome}{$_}) for sort keys %{ $f{genome} };
    }
    for my $name (sort keys %{ $f{vcf} || {} }) {
        write_file(File::Spec->catfile($dir, $name), $f{vcf}{$name});
    }
    for my $rel (sort keys %{ $f{snps_in} || {} }) {
        my $abs = File::Spec->catfile($dir, 'snps_in', $rel);
        make_path(dirname($abs));
        write_file($abs, $f{snps_in}{$rel});
    }
    push @written, $f{name};
}

sub write_file {
    my ($path, $content) = @_;
    open my $fh, '>', $path or die "write $path: $!\n";
    print {$fh} $content;
    close $fh;
}

### --- VCF construction ------------------------------------------------------------------------
###
### rows are [position, [genotype, filter] per sample, alt_offset]. REF comes from $seq; the ALT list
### is built so that sample-N's genotype index selects a base that differs from REF.
sub vcf {
    my (%a) = @_;
    my @contigs = @{ $a{contigs} };
    my @samples = @{ $a{samples} };
    my $seq     = $a{seq};
    my $format  = $a{format} || 'GT:PL:DP:FI';
    my $chrom   = defined $a{chrom} ? $a{chrom} : $contigs[0];

    my $out = "##fileformat=VCFv4.2\n";
    $out .= "##contig=<ID=$_,length=" . length($seq) . ">\n" for @contigs;
    $out .= join("\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                       $a{info_name} || 'INFO', $a{format_name} || 'FORMAT', @samples) . "\n";

    my @asserts;
    for my $row (@{ $a{rows} }) {
        my ($pos, $calls, $alts, $info) = @$row;
        my $ref = substr($seq, $pos - 1, 1);
        push @asserts, [$pos, $ref];
        my @alt = map { alt_for($ref, $_) } @{ $alts || [0] };
        my @fields;
        for my $c (@$calls) {
            my ($gt, $fi) = @$c;
            ### The FORMAT keys are GT:PL:DP:FI, so the sample field has to carry four colon-separated
            ### values in that order.
            push @fields, "$gt:0,10,20:5:$fi";
        }
        $out .= join("\t", $chrom, $pos, '.', $ref, join(',', @alt), '.', 'PASS',
                     defined $info ? $info : 'AC=2', $format, @fields) . "\n";
    }

    ### The assertion the header comment promises: every REF written must equal the reference base.
    for my $a (@asserts) {
        my ($pos, $ref) = @$a;
        my $actual = substr($seq, $pos - 1, 1);
        die "REF mismatch at position $pos: wrote '$ref', reference has '$actual'\n"
            unless $ref eq $actual;
    }
    return $out;
}

### A five-column SNP annotation of the kind --skip_filtering expects, with the >chr header line
### filter_relevant_SNP_calls_from_VCF writes at :930.
sub snp_file {
    my ($chrom, @rows) = @_;
    my $out = ">$chrom\n";
    my $n   = 0;
    for my $r (@rows) {
        my ($pos, $strand, $allele) = @$r;
        $out .= join("\t", ++$n, $chrom, $pos, $strand, $allele) . "\n";
    }
    return $out;
}

###############################################################################
### The fixtures
###############################################################################

my @STD_POS  = (30, 60, 130, 200, 250, 300);
my $STD_ARGS = '--vcf_file snps.vcf --reference_genome genome --strain STRAIN_A';

### Six homozygous high-confidence SNPs in STRAIN_A, none in STRAIN_B. Enough records that a
### hash-ordered output file visibly reorders between runs, which the >=5-record rule requires.
my $STD_VCF = vcf(
    contigs => ['1'],
    samples => ['STRAIN_A', 'STRAIN_B'],
    seq     => $SEQ,
    rows    => [ map { [ $_, [ ['1/1', 1], ['0/0', 1] ] ] } @STD_POS ],
);

my %STD_GENOME = ('1.fa' => fasta('1', $SEQ));

### --- Modes ------------------------------------------------------------------------------------

fixture(
    name   => 'nmask_basic',
    args   => $STD_ARGS,
    readme => 'Default single-strain N-masking: VCF filtering, genome read, SNP application, both reports, the gzipped all-SNP archive',
    genome => \%STD_GENOME,
    vcf    => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name   => 'full_sequence',
    args   => "$STD_ARGS --full_sequence",
    readme => '--full_sequence does not disable N-masking: both output directories are written',
    genome => \%STD_GENOME,
    vcf    => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name   => 'no_nmasking',
    args   => "$STD_ARGS --no_nmasking",
    readme => '--no_nmasking forces --full_sequence and suppresses the N-masked directory entirely',
    genome => \%STD_GENOME,
    vcf    => { 'snps.vcf' => $STD_VCF },
);

### Every classification in read_snp_files, and five records of each so that hash iteration order can
### actually reorder the output. Below five, a run of identical output across two hash seeds proves
### nothing: with one line per file the order cannot vary at all, and with three it repeats by chance
### one time in six. The five normalisation rules these files exercise are only tested at this size.
my @DUAL_ROWS;
{
    my $pos = 10;
    my @plan = (
        ### [ how many, strain-1 call, strain-2 call, alt indices, what it reaches ]
        [ 5, ['1/1', 1], ['1/1', 1], [0],    'same in both -> SNPs_in_common'                     ],
        [ 5, ['1/1', 1], ['2/2', 1], [0, 1], 'both homozygous for a different base -> different'   ],
        [ 5, ['1/1', 1], ['0/0', 1], [0],    'unique to Ref -> STRAIN_A_specific and the prefixed ID' ],
        [ 5, ['0/0', 1], ['1/1', 1], [0],    'unique to SNP -> STRAIN_B_specific'                  ],
        [ 2, ['1/1', 0], ['1/1', 1], [0],    'low confidence in strain 1'                          ],
        [ 2, ['1/1', 1], ['1/1', 0], [0],    'low confidence in strain 2'                          ],
        [ 1, ['0/1', 1], ['1/1', 1], [0],    'heterozygous in strain 1'                            ],
    );
    for my $p (@plan) {
        my ($n, $a, $b, $alts) = @$p;
        for (1 .. $n) {
            push @DUAL_ROWS, [ $pos, [ $a, $b ], $alts ];
            $pos += 10;
        }
    }
}

fixture(
    name   => 'dual_hybrid',
    args   => "$STD_ARGS --strain2 STRAIN_B --genome_build TESTBUILD",
    readme => 'Dual hybrid: strain comparison, new Ref/SNP annotation, the third genome built on strain 1 full sequence, and --genome_build in six filenames. Five records per classification so hash order can reorder every affected output',
    genome => \%STD_GENOME,
    vcf    => {
        'snps.vcf' => vcf(
            contigs => ['1'],
            samples => ['STRAIN_A', 'STRAIN_B'],
            seq     => $SEQ,
            rows    => \@DUAL_ROWS,
        ),
    },
);

### Genome carries a chromosome the SNP directory has no file for, so read_snps warns and returns an
### empty list and create_modified_chromosome clears the array.
fixture(
    name    => 'skip_filtering',
    args    => '--skip_filtering --reference_genome genome --strain STRAIN_A',
    readme  => '--skip_filtering: hardcoded chromosome list, read_snps standalone, the missing-SNP-file warn and the cleared-array branch',
    genome  => { '1.fa' => fasta('1', $SEQ), 'X.fa' => fasta('X', $SEQ2) },
    snps_in => {
        'SNPs_STRAIN_A/chr1.txt' => snp_file('1',
            map { [ $_, 1, substr($SEQ, $_ - 1, 1) . '/' . alt_for(substr($SEQ, $_ - 1, 1)) ] } @STD_POS),
    },
);

### strand -1 makes read_snps:866 complement both alleles, so the file's ref allele has to be the
### complement of the reference base for create_modified_chromosome:626 to accept it.
my %COMPLEMENT = (A => 'T', C => 'G', G => 'C', T => 'A');
fixture(
    name    => 'reverse_strand',
    args    => '--skip_filtering --reference_genome genome --strain STRAIN_A',
    readme  => 'Reverse-strand SNPs are complemented by read_snps; a non-DNA allele is skipped with a warning. Only reachable through --skip_filtering, since VCF filtering always writes strand 1',
    genome  => { '1.fa' => fasta('1', $SEQ) },
    snps_in => {
        'SNPs_STRAIN_A/chr1.txt' => snp_file('1',
            (map {
                my $ref = substr($SEQ, $_ - 1, 1);
                [ $_, -1, $COMPLEMENT{$ref} . '/' . $COMPLEMENT{ alt_for($ref) } ]
            } @STD_POS),
            [ 310, 1, 'R/Y' ],       # not G, A, T or C
        ),
    },
);

### --- VCF parsing ------------------------------------------------------------------------------

### Every branch of the genotype and filter table in one file.
fixture(
    name   => 'genotypes',
    args   => $STD_ARGS,
    readme => 'The whole filtering table: 0/0, 1/1, 2/2 and 3/3 against multi-ALT records, 0/1 and ./. into other, FI=0, a non-ATCG REF, and an undefined ALT',
    genome => \%STD_GENOME,
    vcf    => {
        'snps.vcf' => do {
            my $v = vcf(
                contigs => ['1'],
                samples => ['STRAIN_A', 'STRAIN_B'],
                seq     => $SEQ,
                rows    => [
                    [  30, [ ['0/0', 1], ['0/0', 1] ] ],                 # same as reference
                    [  60, [ ['1/1', 1], ['0/0', 1] ] ],                 # homozygous, high confidence
                    [  90, [ ['2/2', 1], ['0/0', 1] ], [0, 1] ],         # second alternative allele
                    [ 120, [ ['3/3', 1], ['0/0', 1] ], [0, 1, 2] ],      # third alternative allele
                    [ 150, [ ['0/1', 1], ['0/0', 1] ] ],                 # heterozygous
                    [ 180, [ ['./.', 1], ['0/0', 1] ] ],                 # no genotype call
                    [ 210, [ ['1/1', 0], ['0/0', 1] ] ],                 # low confidence
                    [ 240, [ ['1/1', 1], ['0/0', 1] ] ],                 # a second high-confidence SNP
                ],
            );
            ### Two rows the builder cannot express, because both deliberately break the REF/ALT
            ### contract it asserts: a REF that is not a single ATCG base, and an ALT with no clearly
            ### defined base. Appended verbatim.
            $v .= join("\t", 1, 270, '.', 'CA', 'G', '.', 'PASS', 'AC=2', 'GT:PL:DP:FI',
                       '1/1:0,10,20:5:1', '0/0:0,10,20:5:1') . "\n";
            $v .= join("\t", 1, 300, '.', substr($SEQ, 299, 1), '<DEL>', '.', 'PASS', 'AC=2',
                       'GT:PL:DP:FI', '1/1:0,10,20:5:1', '0/0:0,10,20:5:1') . "\n";
            $v;
        },
    },
);

### Named for the v7 release so process_commandline:1435 turns --v7_VCF on by itself.
fixture(
    name    => 'vcf_v7',
    args    => '--vcf_file mgp_REL2005_snps_indels.vcf.gz --reference_genome genome --strain STRAIN_A',
    readme  => 'v7 VCF layout: --v7_VCF set from the filename, INDEL records skipped, and the extra INDEL line in both reports',
    markers => ['gzip_vcf'],
    genome  => \%STD_GENOME,
    vcf     => {
        'mgp_REL2005_snps_indels.vcf' => do {
            my $v = vcf(
                contigs => ['1'],
                samples => ['STRAIN_A', 'STRAIN_B'],
                seq     => $SEQ,
                rows    => [ map { [ $_, [ ['1/1', 1], ['0/0', 1] ] ] } @STD_POS ],
            );
            ### INDEL lives in INFO, which is what the v7 branch tests.
            $v .= join("\t", 1, 310, '.', substr($SEQ, 309, 1), 'GG', '.', 'PASS', 'INDEL;AC=2',
                       'GT:PL:DP:FI', '1/1:0,10,20:5:1', '0/0:0,10,20:5:1') . "\n";
            $v;
        },
    },
);

fixture(
    name    => 'vcf_gzipped',
    args    => '--vcf_file snps.vcf.gz --reference_genome genome --strain STRAIN_A',
    readme  => 'A gzipped VCF under a name that does not trigger v7 detection, so --v7_VCF stays off',
    markers => ['gzip_vcf'],
    genome  => \%STD_GENOME,
    vcf     => { 'snps.vcf' => $STD_VCF },
);

### --- Genome reading ---------------------------------------------------------------------------

### Chromosomes 10 and 2 pin the string sort - 10 sorts before 2 - and chromosome 3 is in the genome
### but not the VCF, so it takes the unmodified-sequence branch.
fixture(
    name   => 'multi_chrom',
    args   => $STD_ARGS,
    readme => 'Two SNP-carrying chromosomes plus one absent from the VCF: string ordering of chromosome names, the no-SNP-information branch, and the detected-chromosome list',
    genome => {
        '10.fa' => fasta('10', $SEQ),
        '2.fa'  => fasta('2',  $SEQ2),
        '3.fa'  => fasta('3',  $SEQ3),
    },
    vcf => {
        'snps.vcf' => do {
            my $head = "##fileformat=VCFv4.2\n##contig=<ID=10,length=" . length($SEQ)
                     . ">\n##contig=<ID=2,length=" . length($SEQ2) . ">\n"
                     . join("\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                                  'FORMAT', 'STRAIN_A', 'STRAIN_B') . "\n";
            my $body = '';
            for my $c (['10', $SEQ, [30, 60, 130]], ['2', $SEQ2, [40, 90, 150]]) {
                my ($chr, $seq, $pos) = @$c;
                for my $p (@$pos) {
                    my $ref = substr($seq, $p - 1, 1);
                    $body .= join("\t", $chr, $p, '.', $ref, alt_for($ref), '.', 'PASS', 'AC=2',
                                  'GT:PL:DP:FI', '1/1:0,10,20:5:1', '0/0:0,10,20:5:1') . "\n";
                }
            }
            $head . $body;
        },
    },
);

fixture(
    name   => 'multifasta',
    args   => $STD_ARGS,
    readme => 'One FastA file holding two chromosomes: the multi-FastA branch of read_genome_into_memory',
    genome => { 'both.fa' => fasta('1', $SEQ) . fasta('2', $SEQ2) },
    vcf    => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name   => 'fasta_extension',
    args   => $STD_ARGS,
    readme => 'A .fasta reference, reached only when the .fa glob finds nothing',
    genome => { '1.fasta' => fasta('1', $SEQ) },
    vcf    => { 'snps.vcf' => $STD_VCF },
);

### Two files carrying the same header reach the second duplicate-name die, at :1341 - the one whose
### message ends in a full stop. The in-loop die at :1318 needs a third occurrence of the name,
### because the check runs before the previous chromosome has been stored.
fixture(
    name   => 'duplicate_chrom_name',
    args   => $STD_ARGS,
    readme => 'Two FastA files with the same chromosome name abort at the end-of-file duplicate check (the message ending in a full stop), not the in-loop one',
    genome => { 'a.fa' => fasta('1', $SEQ), 'b.fa' => fasta('1', $SEQ2) },
    vcf    => { 'snps.vcf' => $STD_VCF },
);

### A header with no sequence is stored as the empty string, which is falsy, so the chromosome-name
### check fires for a chromosome that is present - reporting it as not found and then listing it.
fixture(
    name   => 'empty_chromosome',
    args   => $STD_ARGS,
    readme => 'A header-only FastA entry triggers the Ensembl-versus-UCSC abort for a chromosome that is present, and the diagnostic lists the name it says was not found',
    genome => { '1.fa' => ">1\n", '2.fa' => fasta('2', $SEQ2) },
    vcf    => {
        'snps.vcf' => do {
            my $head = "##fileformat=VCFv4.2\n##contig=<ID=1,length=320>\n##contig=<ID=2,length="
                     . length($SEQ2) . ">\n"
                     . join("\t", '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                                  'FORMAT', 'STRAIN_A', 'STRAIN_B') . "\n";
            my $ref = substr($SEQ2, 39, 1);
            $head . join("\t", 2, 40, '.', $ref, alt_for($ref), '.', 'PASS', 'AC=2', 'GT:PL:DP:FI',
                         '1/1:0,10,20:5:1', '0/0:0,10,20:5:1') . "\n";
        },
    },
);

### --- Errors and edges ------------------------------------------------------------------------

### The VCF names chromosome 1; the reference calls it chr1. Nothing matches at all, which the
### whole-genome check catches before any sequence is written.
fixture(
    name   => 'chrom_name_mismatch',
    args   => $STD_ARGS,
    readme => 'Ensembl-style VCF chromosome names against UCSC-style reference names abort with both name lists, before an unmodified genome can be written',
    genome => { 'chr1.fa' => fasta('chr1', $SEQ) },
    vcf    => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name    => 'poison_gzip',
    args    => $STD_ARGS,
    readme  => 'A failing gzip aborts the run: the close on the pipe is the only place it can be noticed, and SNPsplit reads the file it would have written',
    markers => ['poison_gzip'],
    genome  => \%STD_GENOME,
    vcf     => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name    => 'missing_genome',
    args    => '--vcf_file snps.vcf --strain STRAIN_A',
    readme  => 'A missing --reference_genome aborts with a non-zero exit status',
    markers => ['allow_empty'],
    vcf     => { 'snps.vcf' => $STD_VCF },
);

### A body record on a contig with no ##contig line has no filehandle, so the print dies with a raw
### perl diagnostic rather than a message naming the problem.
fixture(
    name   => 'undeclared_contig',
    args   => $STD_ARGS,
    readme => 'A VCF record on a chromosome with no ##contig header aborts with a raw perl diagnostic instead of a diagnosis',
    genome => \%STD_GENOME,
    vcf    => {
        'snps.vcf' => do {
            my $v = $STD_VCF;
            $v =~ s/^##contig.*\n//m;
            $v;
        },
    },
);

### The column has to be renamed rather than removed: removing it leaves nine columns, detect_strains
### skips indices up to 8, and the run dies earlier on an empty strain list instead.
fixture(
    name   => 'no_format_column',
    args   => $STD_ARGS,
    readme => 'A renamed FORMAT column aborts in the header parser. Renamed rather than deleted, because deleting it leaves nine columns and the run dies earlier on an empty strain list',
    genome => \%STD_GENOME,
    vcf    => {
        'snps.vcf' => vcf(
            contigs     => ['1'],
            samples     => ['STRAIN_A', 'STRAIN_B'],
            seq         => $SEQ,
            format_name => 'FORMATX',
            rows        => [ map { [ $_, [ ['1/1', 1], ['0/0', 1] ] ] } @STD_POS ],
        ),
    },
);

fixture(
    name   => 'no_fi_field',
    args   => $STD_ARGS,
    readme => 'A FORMAT field without FI aborts when the filter index is looked up',
    genome => \%STD_GENOME,
    vcf    => {
        'snps.vcf' => do {
            my $v = vcf(
                contigs => ['1'],
                samples => ['STRAIN_A', 'STRAIN_B'],
                seq     => $SEQ,
                format  => 'GT:PL:DP',
                rows    => [ map { [ $_, [ ['1/1', 1], ['0/0', 1] ] ] } @STD_POS ],
            );
            ### The builder always writes four sample values; drop the fourth to match the FORMAT.
            $v =~ s/:0,10,20:5:[01]\b/:0,10,20:5/g;
            $v;
        },
    },
);

### git cannot store an empty directory, so the placeholder is what makes the fixture exist at all.
### The .fa and .fasta globs both find nothing, which is the point.
fixture(
    name   => 'empty_genome_folder',
    args   => $STD_ARGS,
    readme => 'A reference folder with no FastA files aborts. The .gitkeep placeholder exists because git cannot store an empty directory',
    genome => { '.gitkeep' => '' },
    vcf    => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name    => 'strain_not_in_vcf',
    args    => '--vcf_file snps.vcf --reference_genome genome --strain NOPE',
    readme  => 'An unknown strain name lists the available strains and aborts',
    markers => ['allow_empty'],
    genome  => \%STD_GENOME,
    vcf     => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name    => 'list_strains',
    args    => '--vcf_file snps.vcf --list_strains',
    readme  => '--list_strains prints the sample names to stdout and exits 0',
    markers => ['allow_empty'],
    vcf     => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name    => 'no_snp_folder',
    args    => '--skip_filtering --reference_genome genome --strain STRAIN_A',
    readme  => '--skip_filtering with no SNP folder aborts with the message telling the user to drop the option. The only guidance this failure gives',
    markers => ['allow_empty'],
    genome  => \%STD_GENOME,
);

fixture(
    name   => 'not_fasta',
    args   => $STD_ARGS,
    readme => 'A .fa file whose first line is not a header aborts. The only check on reference file format',
    genome => { '1.fa' => "this is not a FastA header\nACGTACGTAC\n" },
    vcf    => { 'snps.vcf' => $STD_VCF },
);

fixture(
    name    => 'vcf_missing',
    args    => '--vcf_file absent.vcf --reference_genome genome --strain STRAIN_A',
    readme  => 'A named VCF that does not exist aborts before anything else is read',
    markers => ['allow_empty'],
);

fixture(
    name    => 'same_strain',
    args    => '--vcf_file snps.vcf --reference_genome genome --strain STRAIN_A --strain2 STRAIN_A',
    readme  => 'Strain and strain 2 naming the same sample aborts',
    markers => ['allow_empty'],
    genome  => \%STD_GENOME,
    vcf     => { 'snps.vcf' => $STD_VCF },
);

### A dual hybrid built entirely from pre-made SNP files, with no VCF read at all. Reaches
### read_snp_files against committed archives and the third genome built on strain 1's full sequence,
### neither of which any other --skip_filtering fixture touches.
###
### The .gz files here hold plain text; the runner compresses them on the way in.
{
    my @a_pos = (30,  60,  90, 120, 150);
    my @b_pos = (30,  60, 200, 230, 260);
    my $rows  = sub {
        my @r = map { my $r = substr($SEQ, $_ - 1, 1); [ $_, 1, "$r/" . alt_for($r) ] } @{ $_[0] };
        return @r;
    };
    my $plain = sub {
        my ($n, $out) = (0, '');
        for my $r (@{ $_[0] }) {
            $out .= join("\t", ++$n, 1, $r->[0], $r->[1], $r->[2]) . "\n";
        }
        return $out;
    };
    my @a_rows = $rows->(\@a_pos);
    my @b_rows = $rows->(\@b_pos);

    fixture(
        name    => 'skip_filtering_dual',
        args    => '--skip_filtering --reference_genome genome --strain STRAIN_A --strain2 STRAIN_B',
        readme  => '--dual_hybrid under --skip_filtering, built from committed SNP files and archives with no VCF: --strain2 promotes to --dual_hybrid, which sets --full_sequence for the third genome to read back',
        genome  => { '1.fa' => fasta('1', $SEQ) },
        snps_in => {
            'SNPs_STRAIN_A/chr1.txt'            => snp_file('1', @a_rows),
            'SNPs_STRAIN_B/chr1.txt'            => snp_file('1', @b_rows),
            'all_SNPs_STRAIN_A_GRCm39.txt.gz'   => $plain->(\@a_rows),
            'all_SNPs_STRAIN_B_GRCm39.txt.gz'   => $plain->(\@b_rows),
        },
    );
}

### An existing SNP folder suppresses the creating-it-for-you notice, and an existing archive triggers
### the overwrite notice. The only two messages about clobbering the user's own files.
fixture(
    name         => 'rerun_overwrite',
    args         => $STD_ARGS,
    readme       => 'A second run over its own output: the SNP folder already exists so no notice is given, and the all-SNP archive is reported as overwritten',
    genome       => \%STD_GENOME,
    vcf          => { 'snps.vcf' => $STD_VCF },
    pre_existing => ['all_SNPs_STRAIN_A_GRCm39.txt.gz'],
    snps_in      => { 'SNPs_STRAIN_A/chr1.txt' => ">1\n" },
);

printf "Wrote %d fixtures to %s\n", scalar @written, $fixtures;
print "  $_\n" for @written;
print "\nNow record the expected output and read the diff:\n  test/run_genome_tests.pl --update\n";
