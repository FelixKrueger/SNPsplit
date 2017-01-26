## v0.3.1_dev

SNPsplit
-----

- Added a check to `genome_genome_preparation` that produces a [FATAL ERROR] if the stored chromosome names are not the same as the ones in the VCF file (which is a rather common mistake when people use the Ensembl VCF file but get the genome from UCSC. This should change soon if and when Ensembl adopts the same standard used by NCBI/UCSC).

- Changed the `samtools` command throughout SNPsplit to now correctly use the path supplied by the user with `--samtools_path`. Thanks to Kenzo Hillion for spotting this (see here 77286e1e1ad686ef5e6efcdbda826fcd5e4a5ce2). 

- Option `--genome_build [NAME]` should now work as intended (used to be `--build` only).

SNPsplit_genome_preparation
------

- Added a new version of the genome preparation script that can deal with the latest version of the VCF file for the old NCBIM37 genome build ("mgp.v2.snps.annot.reformat.vcf.gz"). The script is called "SNPsplit_genome_preparation_v2VCF" and may be found in the folder "outdated_VCF_versions" on Github.
