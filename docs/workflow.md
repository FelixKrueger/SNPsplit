# The SNPsplit workflow

The general workflow for working with SNPsplit is:

1. Obtain the latest SNP annotation file from the [Mouse Genomes Project](https://www.mousegenomes.org/)

2. Download the Mouse reference genome from Ensembl ([GRCm39](https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/))

3. Create an N-masked reference genome using the VCF file from 1) (either single strain or dual hybrid). This step creates the SNP file required for the `SNPsplit` step

4. Use aligner of choice to index the new N-masked reference genome

5. Run alignments against the N-masked reference using your favourite aligner

6. Run `SNPsplit` on the resulting BAM file, using the SNP file generated in 3)


Once this workflow has completed, proceed to downstream analysis of your choice (not part of SNPsplit).
