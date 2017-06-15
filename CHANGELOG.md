## v0.3.3

SNPsplit
-----

- Changed `FindBin qw($Bin)` to `FindBin qw($RealBin)` so that symlinks to `tag2sort` are resolved properly.

- in certain cases, some specific SNPs were only used for the allele assignment if they were methylated. In more detail: In cases where the SNP was either C/G (REF/ALT) or G/C (REF/ALT), and the read was on the opposing strand, only the methylated form of the C on the reverse strand had previously been allowed as a valid expected base. This has now been changed so that both G and A are considered valid for the strain containing a G at the SNP position (see also this [issue](https://github.com/FelixKrueger/SNPsplit/issues/11)).

- Changed the way in which C>T SNPs are handled in the allele-tagging report (note that this was merely a report/interpretation thing and did not have any effect the on the actual results). Previously, reads without a call for genome 1 or genome 2 had been listed as: 
_reads did not contain one of the expected bases at known SNP positions_. 
In a bisulfite setting this also included C>T SNPs however, and hence the number could have been rather high (>10%). I have now changed this so that reads which had at least one C>T SNP and were unassignable at the same time are scored differently:
_reads that were unassignable contained C>T SNPs preventing the assignment_



## v0.3.2

SNPsplit
-----

- Changed the `samtools` command throughout SNPsplit to now correctly use the path supplied by the user with `--samtools_path`. Thanks to Kenzo Hillion for spotting this (see [here](https://github.com/FelixKrueger/SNPsplit/commit/77286e1e1ad686ef5e6efcdbda826fcd5e4a5ce2)). 

- Option `--genome_build [NAME]` should now work as intended (used to be `--build` only).

SNPsplit_genome_preparation
------

- Relaxed SNP filtering criteria to now support multiple homozygous variants for the same position in the genome. This step should incresae the number of usable SNPs slightly (but noticably). See [here](https://github.com/FelixKrueger/SNPsplit/issues/8)

- Changed the SNP filtering for `--dual_hybrid` mode to only include positions where both strains had a high confidence call (irrespective of the nature of the call). This step should greatly reduce the number of false positive allele calls. See [here](https://github.com/FelixKrueger/SNPsplit/issues/9) for more details.

- Added a check to `SNPsplit_genome_preparation` that produces a [FATAL ERROR] if the stored chromosome names are not the same as the ones in the VCF file (which is a rather common mistake when people use the Ensembl VCF file but get the genome from UCSC. This should change soon if and when Ensembl adopts the same standard used by NCBI/UCSC).

- Added a new version of the genome preparation script that can deal with the latest version of the VCF file for the old NCBIM37 genome build ("mgp.v2.snps.annot.reformat.vcf.gz"). The script is called "SNPsplit_genome_preparation_v2VCF" and may be found in the folder "outdated_VCF_versions" on Github. Please note that this does not include the changes to we made the current version (see above).


## v0.3.1 (18-07-2016)

- Manual: Added a fairly detailed section about how SNPs are filtered and processed during the SNPsplit genome preparation so it can be adapted more easily for different VCF files
