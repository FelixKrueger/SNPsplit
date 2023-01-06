# Note for using a UCSC/NCBI genome in conjunction with the VCF file from the Mouse Genomes Project 

Several users have run into problems when using genomes from UCSC or NCBI in conjunction with the VCF file from the Mouse Genomes Project (MGP, https://www.mousegenomes.org/).
The reason for this is that the MGP uses chromosomal coordinates from Ensembl (i.e. `1, 2, 3, X, MT`) whereas UCSC uses chromosome names that look like this: `chr1, chr2, chr3, chrX, chrM`.

We have recently added a check to the SNPsplit genome preparation script that will bail if a chromosome name discrepancy is detected (https://github.com/FelixKrueger/SNPsplit/issues/4). It is however possible to convert the VCF file into a UCSC compatible version by
(a) changing the chromosome name from e.g. `1` to `chr1` and (b) adding changing the chromosome names in the ID field of the VCF file headers. It is normally not necessary to change the name of the mitochondrium from `MT` to `chrM` because no SNP positions are recorded for the MT anyway.

Here is a one line `awk` script that does an Ensembl=>UCSC conversion, but you could of course also run an equivalent script in Python or Perl...
```
awk '{if($1 ~ "^#") {gsub("contig=<ID=", "contig=<ID=chr"); gsub("contig=<ID=chrMT", "contig=<ID=chrM"); print} else {gsub("^MT", "M"); print "chr"$0}}' mgp_REL2021_snps.vcf.gz
```
