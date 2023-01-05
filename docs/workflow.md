---
hide:
  - navigation
---

# The SNPsplit workflow

The general workflow for working with SNPsplit is outlined below (please note that steps 1-3 are only required as a one-off):

1. Obtain the SNP annotation file from the [Mouse Genomes Project](https://www.mousegenomes.org/). The latest file is [v8 (12/2021)](https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz).

2. Download the Mouse reference genome from Ensembl ([GRCm39](https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/))

3. Create an N-masked reference genome using the VCF file from 1.) (either single strain or dual hybrid). This step creates the SNP file required for the `SNPsplit` step

4. Use aligner of choice to index the new N-masked reference genome

5. Run alignments against the N-masked reference using your favourite aligner

6. Run `SNPsplit` on the resulting BAM file, using the SNP file generated in 3.)


Once this workflow has completed, proceed to downstream analysis of your choice (outside scope of SNPsplit).

If you are still using the older mouse genome build GRCm38 you may find some information [here](./genome_prep/legacy.md).
