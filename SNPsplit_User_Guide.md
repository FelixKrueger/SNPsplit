
 


 
# Examples

## Paired-end report (2x50bp):

```
Input file:					'FVBNJ_Cast.bam'
Writing allele-flagged output file to:	'FVBNJ_Cast.allele_flagged.bam'


Allele-tagging report
=====================
Processed 194564995 read alignments in total
149380724 reads were unassignable (76.78%)
35143075 reads were specific for genome 1 (18.06%)
9860248 reads were specific for genome 2 (5.07%)
118662 reads did not contain one of the expected bases at known SNP positions (0.06%)
180948 contained conflicting allele-specific SNPs (0.09%)


SNP coverage report
===================
N-containing reads:	45276050
non-N:			149262062
total:			194564995
Reads had a deletion of the N-masked position (and were thus dropped):	26883 (0.01%)
Of which had multiple deletions of N-masked positions within the same read:	30

Of valid N containing reads,
N was present in the list of known SNPs:		61087551 (99.99%)
N was not present in the list of SNPs:		4773 (0.01%)

Input file:						'FVBNJ_Cast.allele_flagged.bam'
Writing unassigned reads to:				'FVBNJ_Cast.unassigned.bam'
Writing genome 1-specific reads to:			'FVBNJ_Cast.genome1.bam'
Writing genome 2-specific reads to:			'FVBNJ_Cast.genome2.bam'


Allele-specific paired-end sorting report
=========================================
Read pairs/singletons processed in total:		98215744
	thereof were read pairs:			96349251
	thereof were singletons:			1866493
Reads were unassignable (not overlapping SNPs):	61174812 (62.29%)
	thereof were read pairs:	59662537
	thereof were singletons:	1512275
Reads were specific for genome 1:			28657857 (29.18%)
	thereof were read pairs:	28446094
	thereof were singletons:	211763
Reads were specific for genome 2:			8122687 (8.27%)
	thereof were read pairs:	7985424
	thereof were singletons:	137263
Reads contained conflicting SNP information:		260388 (0.27%)
	thereof were read pairs:	255196
	thereof were singletons:	5192
```
 
## Hi-C report (2x100bp):

```
Input file:					Black6_129S1.bam
Writing allele-flagged output file to:	Black6_129S1.allele_flagged.bam


Allele-tagging report
=====================
Processed 94887256 read alignments in total
59662038 reads were unassignable (62.88%)
19851697 reads were specific for genome 1 (20.92%)
15047281 reads were specific for genome 2 (15.86%)
47261 reads did not contain one of the expected bases at known SNP positions (0.05%)
326240 contained conflicting allele-specific SNPs (0.34%)


SNP coverage report
===================
N-containing reads:	35231977
non-N:			59614777
total:			94887256
Reads had a deletion of the N-masked position (and were thus dropped):	40502 (0.04%)
Of which had multiple deletions of N-masked positions within the same read:	59

Of valid N containing reads,
N was present in the list of known SNPs:	57101748 (99.99%)
N was not present in the list of SNPs:	4211 (0.01%)

Input file:						Black6_129S1.allele_flagged.bam'
Writing unassigned reads to:				Black6_129S1.UA_UA.bam'
Writing genome 1-specific reads to:			Black6_129S1.G1_G1.bam'
Writing genome 2-specific reads to:			Black6_129S1.G2_G2.bam'
Writing G1/UA reads to:				Black6_129S1.G1_UA.bam'
Writing G2/UA reads to:				Black6_129S1.G2_UA.bam'
Writing G1/G2 reads to:				Black6_129S1.G1_G2.bam'


Allele-specific paired-end sorting report
=========================================
Read pairs processed in total:			47443628
Read pairs were unassignable (UA/UA):			18862725 (39.76%)
Read pairs were specific for genome 1 (G1/G1):	3533932 (7.45%)
Read pairs were specific for genome 2 (G2/G2):	2592040 (5.46%)
Read pairs were a mix of G1 and UA:			12306421 (25.94%). Of these,
			were G1/UA: 6018598
			were UA/G1: 6287823
Read pairs were a mix of G2 and UA:			9430675 (19.88%). Of these,
			were G2/UA: 4603429
			were UA/G2: 4827246
Read pairs were a mix of G1 and G2:			395296 (0.83%). Of these,
			were G1/G2: 198330
			were G2/G1: 196966
Read pairs contained conflicting SNP information:	322539 (0.68%)
```


 
## BS-Seq report (2x100bp):

```
Input file:					'129_Cast_bismark_bt2_pe.bam'
Writing allele-flagged output file to:	'129_Cast_bismark_bt2_pe.allele_flagged.bam'


Allele-tagging report
=====================
Processed 162441396 read alignments in total
Reads were unaligned and hence skipped: 0 (0.00%)
109109113 reads were unassignable (67.17%)
30267901 reads were specific for genome 1 (18.63%)
22697499 reads were specific for genome 2 (13.97%)
15807753 reads did not contain one of the expected bases at known SNP positions (9.73%)
366883 contained conflicting allele-specific SNPs (0.23%)


SNP coverage report
===================
SNP annotation file:	../all_Cast_SNPs_129S1_reference.mgp.v4.txt.gz
N-containing reads:	68984287
non-N:			93301360
total:			162441396
Reads had a deletion of the N-masked position (and were thus dropped):	155749 (0.10%)
Of which had multiple deletions of N-masked positions within the same read:	65

Of valid N containing reads,
N was present in the list of known SNPs:	119119643 (99.99%)
Positions were skipped since they involved C>T SNPs:	38464451
N was not present in the list of SNPs:		7517 (0.01%)

Input file:						129_Cast_bismark_bt2_pe.allele_flagged.bam'
Writing unassigned reads to:				129_Cast_bismark_bt2_pe.unassigned.bam'
Writing genome 1-specific reads to:			129_Cast_bismark_bt2_pe.genome1.bam'
Writing genome 2-specific reads to:			129_Cast_bismark_bt2_pe.genome2.bam'


Allele-specific paired-end sorting report
=========================================
Read pairs/singletons processed in total:		81220698
	thereof were read pairs:			81220698
	thereof were singletons:			0
Reads were unassignable (not overlapping SNPs):		40420625 (49.77%)
	thereof were read pairs:	40420625
	thereof were singletons:	0
Reads were specific for genome 1:			23037433 (28.36%)
	thereof were read pairs:	23037433
	thereof were singletons:	0
Reads were specific for genome 2:			17303663 (21.30%)
	thereof were read pairs:	17303663
	thereof were singletons:	0
Reads contained conflicting SNP information:		458977 (0.57%)
	thereof were read pairs:	458977
	thereof were singletons:	0
```
 
## Full list of options for SNPsplit


**USAGE:** `SNPsplit [options] --snp_file <SNP.file.gz> [input file(s)]`



`Input file(s)` 
- Mapping output file in SAM or BAM format. SAM files (ending in `.sam`) will first be converted to BAM files.

`--snp_file`
- Mandatory file specifying SNP positions to be considered, may be a plain text file of gzip compressed. Currently, the SNP file is expected to be in the following format:
```
   SNP-ID     Chromosome  Position    Strand   Ref/SNP
 33941939           9             68878541       1           T/G
```

	Only the information contained in fields 'Chromosome', 'Position' and 'Ref/SNP base' are being used for analysis. The genome referred to as 'Ref' will be used as genome 1, the genome containing the 'SNP' base as genome 2.

`--paired`
- Paired-end mode. (Default: OFF).

`-o/--outdir <dir>`
- Write all output files into this directory. By default the output files will be written into the same folder as the input file(s). If the specified folder does not exist, SNPsplit will attempt to create it first. The path to the output folder can be either relative or absolute.


`--singletons`
- If the allele-tagged paired-end file also contains singleton alignments (which is the default for e.g. TopHat), these will be written out to extra files (ending in `_st.bam`) instead of writing everything to combined paired-end and singleton files. Default: OFF.

`--no_sort`
- This option skips the sorting step if BAM files are already sorted by read name (e.g. Hi-C files generated by [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/)). Please note that setting `--no_sort` for unsorted paired-end  files will break the tagging process!

`--hic`
- Assumes Hi-C data processed with [HiCUP](www.bioinformatics.babraham.ac.uk/projects/hicup/) as input, i.e. the input BAM file is paired-end and Reads 1 and 2 follow each other. Thus, this option also sets the flags `--paired` and `--no_sort`. Default: OFF.

`--bisulfite`
- Assumes Bisulfite-Seq data processed with [Bismark](https://github.com/FelixKrueger/Bismark) as input. In paired-end mode (`--paired`), Read 1 and Read 2 of a pair are expected to follow each other in consecutive lines. SNPsplit will run a quick check at the start of a run to see if the provided file appears to be a Bismark file, and set the flags `--bisulfite` and/or `--paired` automatically. In addition it will perform a quick check to see if a paired-end file appears to have been positionally sorted, and if not will set the flag `--no_sort`.

`--samtools-path`
- The path to your Samtools installation, e.g. `/home/user/samtools/`. Does not need to be specified explicitly if Samtools is in the `PATH` environment already.

## SNPsplit-sort specific options (tag2sort):

`--sam`
- The output will be written out in SAM format instead of the default BAM format. SNPsplit will attempt to use the path to Samtools that was specified with `--samtools_path`, or, if it hasn't been specified, attempt to find Samtools in the `PATH` environment. 

`--conflicting/--weird`
- Reads or read pairs that were classified as 'Conflicting' (`XX:Z:CF`) will be written to an extra file (ending in `.conflicting.bam`) instead of being simply skipped. Reads may be classified as 'Conflicting' if a single read contains SNP information for both genomes at the same time, or if the SNP position was deleted from the read. Read-pairs are considered 'Conflicting' if either read is was tagged with the `XX:Z:CF` flag. Default: OFF.

`--help`
- Displays this help information and exits.

`--verbose`
(Very!) verbose output (for debugging).

`--version`
- Displays version information and exits.
 
## Full list of options for SNPsplit_genome_preparation


**USAGE:** `SNPsplit_genome_preparation [options] --vcf_file <file> --reference_genome /path/to/genome/ --strain <strain_name>`


`--vcf_file <file>`
- Mandatory file specifying SNP information for mouse strains from the [Mouse Genomes Project](http://www.sanger.ac.uk/science/data/mouse-genomes-project). The file approved and tested is called `mgp.v5.merged.snps_all.dbSNP142.vcf.gz`. Please note that future versions of this file or entirely different VCF files might not work out-of-the-box but may require some tweaking. SNP calls are read from the VCF files, and high confidence SNPs are written into a folder in the current working directory called `SNPs_<strain_name>/chr<chromosome>.txt`, in the following format:
```
        SNP-ID     Chromosome  Position    Strand   Ref/SNP
example:   33941939        9       68878541       1       T/G
```

`--strain <strain_name>`
- The strain you would like to use as SNP (ALT) genome. Mandatory. For an overview of strain names just run SNPsplit_genome_preparation selecting `--list_strains`.

`--list_strains`
- Displays a list of strain names present in the VCF file for use with `--strain <strain_name>`.

`--dual_hybrid`
- Optional. The resulting genome will no longer relate to the original reference specified with `--reference_genome`. Instead the new Reference (Ref) is defined by `--strain <name>` and the new SNP genome is defined by `--strain2 <strain_name>`. `--dual_hybrid` automatically sets `--full_sequence`.

	**This will invoke a multi-step process:**
**1)** Read/filter SNPs for first strain (specified with `--strain < name>`)
**2)** Write full SNP incorporated (and optionally N-masked) genome sequence for first strain
**3)** Read/filter SNPs for second strain (specified with `--strain2 <name>`)
**4)** Write full SNP incorporated (and optionally N-masked) genome sequence for second strain
**5)** Generate new Ref/Alt SNP annotations for Strain1/Strain2
**6)** Set first strain as new reference genome and construct full SNP incorporated (and optionally N-masked) genome sequences for Strain1/Strain2

`--strain2 <strain_name>`
- Optional for constructing dual hybrid genomes (see `--dual_hybrid` for more information). For an overview of strain names just run SNPsplit_genome_preparation selecting `--list_strains`.

`--reference_genome`
- The path to the reference genome, typically the strain 'Black6' (C57BL/6J), e.g. `--reference_genome /bi/scratch/Genomes/Mouse/GRCm38/`. Expects one or more FastA files in this folder (file extension: `.fa` or `.fasta`).

`--skip_filtering`
- This option skips reading and filtering the VCF file. This assumes that a folder named `SNPs_<Strain_Name>` exists in the working directory, and that text files with SNP information are contained therein in the following format:
```
              SNP-ID     Chromosome  Position    Strand   Ref/SNP
example:   33941939        9       68878541       1       T/G
```

`--nmasking`
- Write out a genome version for the strain specified where Ref bases are replaced with `N`. In the Ref/SNP example T/G the N-masked genome would now carry an N instead of the T. The N-masked genome is written to a folder called  `<strain_name>_N-masked/`. Default: ON.

`--full_sequence`
- Write out a genome version for the strain specified where Ref bases are replaced with the SNP base. In the Ref/SNP example T/G the full sequence genome would now carry a G instead of the T. The full sequence genome is written out to folder called `<strain_name>_full_sequence/`. May be set in addition to `--nmasking`. Default: OFF.

`--no_nmasking`
- Disable N-masking if it is not desirable. Will automatically set `--full_sequence` instead.

`--genome_build <name>`
- Name of the mouse genome build, e.g. mm10. Will be incorporated into some of the output files. Defaults to `GRCm38`.

`--help`
- Displays this help information and exits.

`--version`
- Displays version information and exits.
