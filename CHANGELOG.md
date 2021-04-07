## v0.4.0dev (07 04 2021)

#### SNPsplit_genome_preparation

- Added option `--v7_MGP`; now also accepts the v7 file (`mgp_REL2005_snps_indels.vcf.gz`) of Mouse Genomes Project which may be downloaded here: ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz. INDEL variants are being skipped (this is noted in the report). This new version adds a number of additional strains to choose from, now amounting to 51 strains in total:

**Available genomes to choose from are:**

```LG_J
SEA_GnJ
SM_J
ST_bJ
CAST_EiJ
BALB_cByJ
NON_LtJ
FVB_NJ
RIIIS_J
CE_J
NZO_HlLtJ
C58_J
BTBR_T+_Itpr3tf_J
MOLF_EiJ
BUB_BnJ
C57L_J
CZECHII_EiJ
C57BL_10J
B10.RIII
AKR_J
C3H_HeJ
LP_J
DBA_2J
QSi3
ZALENDE_EiJ
A_J
PL_J
129S1_SvImJ
NZW_LacJ
PWK_PhJ
C57BL_10SnJ
C57BR_cdJ
QSi5
C57BL_6NJ
SWR_J
MA_MyJ
C3H_HeH
SPRET_EiJ
LEWES_EiJ
WSB_EiJ
129P2_OlaHsd
CBA_J
SJL_J
BALB_cJ
KK_HiJ
JF1_MsJ
NZB_B1NJ
I_LnJ
DBA_1J
129S5SvEvBrd
NOD_ShiLtJ
RF_J
```
If the file `mgp_REL2005_snps_indels.vcf.gz` is given, `--v7` is set automatically.

- now attempts to extract the fields `FORMAT` and `INFO` from the VCF file automatically, to get access to the required information `GT` (genotype) and `FI` (filter). See more here: #47.


## v0.4.0 (29 09 2020)

- SNPsplit now supports soft-clipping of reads (`CIGAR` operation `S`). 

- SNPsplit now writes important statistics out in YAML format to enable easier integration into `MultiQC`. If `tag2sort` is called via `SNPsplit` itself, the `...sort.yaml` file will be integrated into the main `...SNPsplit_report.yaml` file (and deleted afterwards)

- Added option `--skip_tag2sort` to allow the separation of the allele-tagging and allele-sorting (`tag2sort`) processes. This might be desired to add a de-duplication step such as `markduplicates` or `deduplicate_bismark` for Nextflow pipelines

- For genomes that consist of chromosomes for which SNPs are recorded, and scaffolds for which there are no SNPs, now all chromosomes and scaffolds are printed to both the N-masked and full sequence genomes (see [here](https://github.com/FelixKrueger/SNPsplit/issues/38)).

- added auto-detection of single-end or paired-end files. This avoids accidentally processing paired-end files in single-end mode [see here](https://github.com/FelixKrueger/SNPsplit/issues/27).

- now making use of variable genome_build instead of using GRCm38 invariably


## v0.3.4 (25 May 2018)

 - made the installation changes for conda
 
 - fixed --dir name in tag2sort

- Fixed output-path handling for paired-end and Hi-C mode (was only working for single-end files).

- Added option `-o/--output_dir` to specify an output directory.

## v0.3.3 (15 June 2017)

#### SNPsplit
-----

- Changed `FindBin qw($Bin)` to `FindBin qw($RealBin)` so that symlinks to `tag2sort` are resolved properly.

- In certain cases, specific SNPs were only used for the allele assignment if they were methylated. In more detail: In cases where the SNP was either C/G (REF/ALT) or G/C (REF/ALT), and the read was on the opposing strand, only the methylated form of the C on the reverse strand had previously been allowed as a valid expected base. This has now been changed so that both G and A are considered valid for the strain containing a G at the SNP position (see also this [issue](https://github.com/FelixKrueger/SNPsplit/issues/11)).

- Changed the way in which C>T SNPs are handled in the allele-tagging report (note that this was merely a report/interpretation thing and did not have any effect the on the actual results). Previously, reads without a call for genome 1 or genome 2 had been listed as:
_reads did not contain one of the expected bases at known SNP positions_.
In a bisulfite setting this also included C>T SNPs however, and hence the number could have been rather high (>10%). I have now changed this so that reads which had at least one C>T SNP and were unassignable at the same time are scored differently:
_reads that were unassignable contained C>T SNPs preventing the assignment_

- Changed all instances of `zcat` to `gunzip -c` in `SNPsplit` and `SNPsplit_genome_preparation` to prevent errors on certain OSX platforms


## v0.3.2 (29 March 2017)

#### SNPsplit
-----

- Changed the `samtools` command throughout SNPsplit to now correctly use the path supplied by the user with `--samtools_path`. Thanks to Kenzo Hillion for spotting this (see [here](https://github.com/FelixKrueger/SNPsplit/commit/77286e1e1ad686ef5e6efcdbda826fcd5e4a5ce2)).

- Option `--genome_build [NAME]` should now work as intended (used to be `--build` only).

#### SNPsplit_genome_preparation
------

- Relaxed SNP filtering criteria to now support multiple homozygous variants for the same position in the genome. This step should incresae the number of usable SNPs slightly (but noticably). See [here](https://github.com/FelixKrueger/SNPsplit/issues/8)

- Changed the SNP filtering for `--dual_hybrid` mode to only include positions where both strains had a high confidence call (irrespective of the nature of the call). This step should greatly reduce the number of false positive allele calls. See [here](https://github.com/FelixKrueger/SNPsplit/issues/9) for more details.

- Added a check to `SNPsplit_genome_preparation` that produces a [FATAL ERROR] if the stored chromosome names are not the same as the ones in the VCF file (which is a rather common mistake when people use the Ensembl VCF file but get the genome from UCSC. This should change soon if and when Ensembl adopts the same standard used by NCBI/UCSC).

- Added a new version of the genome preparation script that can deal with the latest version of the VCF file for the old NCBIM37 genome build ("mgp.v2.snps.annot.reformat.vcf.gz"). The script is called "SNPsplit_genome_preparation_v2VCF" and may be found in the folder "outdated_VCF_versions" on Github. Please note that this does not include the changes to we made the current version (see above).


## v0.3.1 (18 July 2016)

- Manual: Added a fairly detailed section about how SNPs are filtered and processed during the SNPsplit genome preparation so it can be adapted more easily for different VCF files

## v0.3.1_dev (03 Nov 2016)

#### SNPsplit
-----

- Minor fixed for path to Samtools when checking files for BS.


#### SNPsplit_genome_preparation
-----

- Option `--genome_build [NAME]` should now work as intended (used to be `--build` only).

- Added a new version of the genome preparation script that can deal with the latest version of the VCF file for the old NCBIM37 genome build (`mgp.v2.snps.annot.reformat.vcf.gz`). The script is called `SNPsplit_genome_preparation_v2VCF` and may be found in the folder `outdated_VCF_versions` on Github.


## v0.3.1 (18 Jul 2016)

- Added a new section to the SNPsplit manual explaining in more detail how the SNP filtering step works and which lines and fields are relevant.



## v0.3.0 (18 May 2016)

#### SNPsplit
-----

- Added the # of SNPs used for the allele-discrimination to the report file to make it easier to spot errors.

- Now removing `CR` and `LF` line endings when reading in the SNP file. For SNP annotation files copied from a Windows machine we saw problems with no allele-specific reads for genome 2 at all which was due to the invisible `\r` character for the SNP call.

- Changed sorting command for BAM files to also work with Samtools versions 1.3+.

- The sorting report for single-end files is now also written to the report files.


#### SNPsplit_genome_preparation
-----

- Added whole new functionality to construct single- or dual-hybrid genomes starting from VCF files which are obtainable from the Mouse Genomes Project (http://www.sanger.ac.uk/science/data/mouse-genomes-project), here is a brief description of what it does:

SNPsplit_genome_preparation is designed to read in a variant call files from the Mouse Genomes Project (e.g. this latest file: ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz) and generate new genome versions where the strain SNPs are either incorporated into the new genome (full sequence) or masked by the ambiguity nucleobase 'N' (**N-masking**).

SNPsplit_genome_preparation may be run in two different modes:

##### Single strain mode:
- The VCF file is read and filtered for high-confidence SNPs in the strain specified with strain
- The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated either as N-masking (default) or full sequence (option `--full_sequence`)

##### Dual strain mode:
- The VCF file is read and filtered for high-confidence SNPs in the strain specified with `--strain <name>`
- The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated as full sequence and optionally as N-masking
- The VCF file is read one more time and filtered for high-confidence SNPs in strain 2 specified with `--strain2 <name>`
- The filtered high-confidence SNP positions of strain 2 are incorporated as full sequence and optionally as N-masking
- The SNP information of strain and strain 2 relative to the reference genome build are compared, and a new Ref/SNP annotation is constructed whereby the new Ref/SNP information will be Strain/Strain2 (and no longer the standard reference genome strain Black6 (C57BL/6J)).The full genome sequence given with `--strain <name>` is read into memory, and the high-confidence SNP positions between Strain and Strain2 are incorporated as full sequence and optionally as N-masking
- The resulting `.fa` files are ready to be indexed with your favourite aligner. Proved and tested aligners include Bowtie2, Tophat, STAR, HISAT2, HiCUP and Bismark. Please note that STAR and HISAT2 may require you to disable soft-clipping, please see the SNPsplit manual more details
- Both the SNP filtering and the genome preparation write out little report files for record keeping



## v0.2.0 (19 August 2015)

- The name of the SNP annotation file is now displayed on screen and written to the report files.

- Added new option `--bisulfite`. This assumes Bisulfite-Seq data processed with Bismark (www.bioinformatics.babraham.ac.uk/projects/bismark/) as input. In paired-end mode (`--paired`), Read 1 and Read 2 of a pair are expected to follow each other in consecutive lines. SNPsplit will run a quick check at the start of a run to see if the file provided file appears to be a Bismark file, and set the flags `--bisulfite` and/or `--paired` automatically. In addition it will perform a quick check to see if a paired-end file appears to have been positionally sorted, and if not will set the `--no_sort` flag.

- Reads having the unmapped FLAG set in the BAM/SAM file (0x4 bit) are now skipped and excluded from the tagging and sorting process.

- Improved file renaming settings when input file was in SAM format (i.e. no longer deletes the input files..). Also changed renaming settings to only change .bam at the end of reads.




## v0.1.3 (04 November 2014)

- SNPsplit v0.1.3 is an initial beta release and as such is still a work in progress. All basic functions are working.


Feedback or bugs can be sent to: felix.krueger@babraham.ac.uk
