---
hide:
  - navigation
---

05 January, 2023

# SNPsplit [v0.6.0dev]
______

## General

SNPsplit is an allele-specific alignment sorter which is designed to read alignment files in SAM/ BAM format and determine the allelic origin of reads that cover known SNP positions. For this to work a library must have been aligned to a genome which had all SNP positions masked by the ambiguity base ’N’, and aligned using aligners that are capable of using a reference genome which contains ambiguous nucleobases. Examples of supported alignment programs are Bowtie2, HISAT2, Bismark, HiCUP, or STAR (for some tips using HISAT2 or STAR alignments please see below). In addition, a list of all known SNP positions between the two different genomes must be provided using the option `--snp_file`. SNP information to generate N-masked genomes needs to be acquired elsewhere, e.g. for different strains of mice you can find variant call files at the [Mouse Genomes Project page](https://www.mousegenomes.org/). 

SNPsplit offers a separate genome preparation step that allows you to generate N-masked (or fully incorporated) SNP genomes for single or dual hybrid strains for all strains of the Mouse Genomes Project. Please see below for further details.

SNPsplit operates in two stages:

 **I) SNPsplit-tag:** 
SNPsplit analyses reads (single-end mode) or read pairs (paired-end mode) for overlaps with known SNP positions, and writes out a tagged BAM file in the same order as the original file. Unsorted paired-end files are sorted by name first.

 **II) SNPsplit-sort:** 
The tagged BAM file is read in again and sorted into allele-specific files. This process may also be run as a stand-alone module on tagged BAM files (tag2sort).

The SNPsplit-tag module determines whether a read can be assigned to a certain allele and appends an additional optional field `XX:Z:` to each read. The tag can be one of the following:

    XX:Z:UA - Unassigned
    XX:Z:G1 - Genome 1-specific
    XX:Z:G2 - Genome 2-specific
    XX:Z:CF - Conflicting

The SNPsplit-sort module `tag2sort` reads in the tagged BAM file and sorts the reads (or read pairs) according to their `XX:Z:` tag (or the combination of tags for paired-end or Hi-C reads) into sub-files.