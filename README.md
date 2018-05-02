
<p align="center"> <img title="SNPsplit" id="logo_img" src="Images/SNPsplit.png" width=300></p>

## Allele-specific alignment sorting


### Note for using a UCSC/NCBI genome in conjunction with the VCF file from the Mouse Genomes Project 

Several users have run into problems when using genomes from UCSC or NCBI in conjunction with the VCF file from the Mouse Genomes Project (MGP, http://www.sanger.ac.uk/science/data/mouse-genomes-project).
The reason for this is that the MGP uses chromosomal coordinates from Ensembl (i.e. `1, 2, 3, X, MT`) whereas UCSC uses chromosome names that look like this: `chr1, chr2, chr3, chrX, chrM`.

We have recently added a check to the SNPsplit genome preparation script that will bail if a chromosome name discrepancy is detected (https://github.com/FelixKrueger/SNPsplit/issues/4). It is however possible to convert the VCF file into a UCSC compatible version by
(a) changing the chromosome name from e.g. `1` to `chr1` and (b) adding changing the chromosome names in the ID field of the VCF file headers. It is normally not necessary to change the name of the mitochondrium from `MT` to `chrM` because no SNP positions are recorded for the MT anyway.

Here is a one line `awk` script that does an Ensembl=>UCSC conversion, but you could of course also run an equivalent script in Perl or Python...
```
awk '{if($1 ~ "^#") {gsub("contig=<ID=", "contig=<ID=chr"); gsub("contig=<ID=chrMT", "contig=<ID=chrM"); print} else {gsub("^MT", "M"); print "chr"$0}}' mgp.v5.merged.snps_all.dbSNP142.vcf
```

## Installation

SNPsplit is written in Perl and is executed from the command line. To install SNPsplit simply download the latest release of the code from the [Releases page](https://github.com/FelixKrueger/SNPsplit/releases) and extract the files into a SNPsplit installation folder.

SNPsplit requires the following tools installed and ideally available in the `PATH`:
- [Samtools](http://samtools.sourceforge.net/)

## Documentation
The SNPsplit documentation can be found here: [SNPsplit User Guide](./SNPsplit_User_Guide.md)

## Links
- SNPsplit publication at F1000 Research:
  * https://f1000research.com/articles/5-1479/v2
  
- Here is a link to the [SNPsplit project site](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/) at the Babraham Institute.

## Credits

SNPsplit was written by Felix Krueger, part of the [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk) group.

<p align="center"> <img title="Babraham Bioinformatics" id="logo_img" src="Images/bioinformatics_logo.png" width=300></p>
