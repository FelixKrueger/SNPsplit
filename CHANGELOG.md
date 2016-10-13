## v0.3.1_dev

- Added a check to `genome_genome_preparation` that produces a [FATAL ERROR] if the stored chromosome names are not the same as the ones in the VCF file (which is a rather common mistake when people use the Ensembl VCF file but get the genome from UCSC. This should change soon if and when Ensembl adopts the same standard used by NCBI/UCSC).

- Changed the `samtools` command throughout SNPsplit to now correctly use the path supplied by the user with `--samtools_path`. Thanks to Kenzo Hillion for spotting this (see here 77286e1e1ad686ef5e6efcdbda826fcd5e4a5ce2). 
