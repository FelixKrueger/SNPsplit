SNPsplit (Rust)

One binary, snpsplit, and three links to it under the classic names:

  SNPsplit                     tag and sort alignments by parental genome
  tag2sort                     sort an already tagged BAM file
  SNPsplit_genome_preparation  build N-masked or SNP-introduced genomes

Put this folder on your PATH, or copy all four (keeping the links) into a folder that is.
It needs nothing else installed: no Perl and no samtools.

On macOS, a downloaded binary may be blocked until you run:
  xattr -d com.apple.quarantine snpsplit

Documentation: https://felixkrueger.github.io/SNPsplit
Licence: LICENSE (GPL-3.0); third-party licences: THIRD-PARTY-NOTICES.md
