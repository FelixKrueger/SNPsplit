# Genome preparation

!!! note
    USAGE:    SNPsplit_genome_preparation  [options] --vcf_file <file> --reference_genome /path/to/genome/ --strain <strain_name>

`--vcf_file <file>`

Mandatory file specifying SNP information for mouse strains from the [Mouse Genomes Project](https://www.mousegenomes.org/): [v8 annotations](https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/). The file used and approved is called 'mgp_REL2021_snps.vcf.gz'. Please note that future versions of this file or entirely different VCF files might not work out-of-the-box but may require some tweaking. SNP calls are read from the VCF files, and high confidence SNPs are written into  a folder in the current working directory called SNPs_<strain_name>/chr<chromosome>.txt, in the following format:

```             SNP-ID     Chromosome  Position    Strand   Ref/SNP
 example:   33941939        9       68878541       1       T/G
```

`--v7_VCF`

This will use the file 'mgp_REL2005_snps_indels.vcf.gz' instead of the mgp v8 file mentioned above, for backward compatibility reasons. This file contains both SNP and INDEL information, but INDELs are skipped. 
!!! note
    NOTE: The v5 and v7 files work for the GRCm38 (now outdated) genome build!

`--strain <strain_name>`

The strain you would like to use as SNP (ALT) genome. Mandatory. For an overview of strain names just run SNPsplit_genome_preparation selecting `--list_strains`

`--list_strains`

- Displays a list of strain names present in the VCF file for use with `--strain <strain_name>`

--dual_hybrid`              
Optional. The resulting genome will no longer relate to the original reference specified with '--reference_genome'.
Instead the new Reference (Ref) is defined by '--strain <strain_name>' and the new SNP genome
is defined by '--strain2 <strain_name>'. '--dual_hybrid' automatically sets '--full_sequence'.
This will invoke a multi-step process:
   1) Read/filter SNPs for first strain (specified with '--strain <strain_name>')
   2) Write full SNP incorporated (and optionally N-masked) genome sequence for first strain
   3) Read/filter SNPs for second strain (specified with '--strain2 <strain_name>')
   4) Write full SNP incorporated (and optionally N-masked) genome sequence for second strain
   5) Generate new Ref/Alt SNP annotations for Strain1/Strain2
   6) Set first strain as new reference genome and construct full SNP incorporated (and optionally
      N-masked) genome sequences for Strain1/Strain2


--strain2 <strain_name>       Optional for constructing dual hybrid genomes (see '--dual_hybrid' for more information). For an
                              overview of strain names just run SNPsplit_genome_preparation selecting '--list_strains'.

--reference_genome            The path to the reference genome, typically the strain 'Black6' (C57BL/6J), e.g.
                              '--reference_genome /scratch/Genomes/Mouse/GRCm39/'. Expects one or more FastA files in this folder
                              (file extension: .fa or .fasta).

--skip_filtering              This option skips reading and filtering the VCF file. This assumes that a folder named
                              'SNPs_<Strain_Name>' exists in the working directory, and that text files with SNP information
                              are contained therein in the following format:

                                          SNP-ID     Chromosome  Position    Strand   Ref/SNP
                              example:   33941939        9       68878541       1       T/G

--nmasking                    Write out a genome version for the strain specified where Ref bases are replaced with 'N'. In the
                              Ref/SNP example T/G the N-masked genome would now carry an N instead of the T. The N-masked genome
                              is written to a folder called  '<strain_name>_N-masked/'. Default: ON.

--full_sequence               Write out a genome version for the strain specified where Ref bases are replaced with the SNP base.
                              In the Ref/SNP example T/G the full sequence genome would now carry a G instead of the T. The full
                              sequence genome is written out to folder called '<strain_name>_full_sequence/'. May be set in
                              addition to '--nmasking'. Default: OFF.

--no_nmasking                 Disable N-masking if it is not desirable. Will automatically set '--full_sequence' instead.

--genome_build [name]         Name of the genome build incorporated into some of the output files. Defaults to 'GRCm39'.

--help                        Displays this help information and exits.

--version                     Displays version information and exits.