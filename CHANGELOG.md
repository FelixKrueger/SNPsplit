## v0.3.2_dev

SNPsplit
-----

- Changed `FindBin qw($Bin)` to `FindBin qw($RealBin)` so that symlinks to `tag2sort` are resolved properly.

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
