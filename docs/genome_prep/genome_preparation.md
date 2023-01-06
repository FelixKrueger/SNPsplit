# SNPsplit genome preparation

`SNPsplit_genome_preparation` is designed to read in a variant call file from the Mouse Genomes Project (this [latest v8 file](https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz)) and generate new genome versions where the strain SNPs are either incorporated into the new genome (full sequence) or masked by the ambiguity nucleobase `N` (**N-masking**).

`SNPsplit_genome_preparation` may be run in two different modes:

### Single strain mode:

**1)** The VCF file is read and filtered for high-confidence SNPs in the strain specified with `--strain <name>`

**2)** The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated either as N-masking (default), or full sequence (option `--full_sequence`)

### Dual strain mode:
**1)** The VCF file is read and filtered for high-confidence SNPs in the strain specified with `--strain <name>`

**2)** The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated as full sequence and optionally as N-masking

**3)** The VCF file is read one more time and filtered for high-confidence SNPs in strain 2 specified with `--strain2 <name>`

**4)** The filtered high-confidence SNP positions of strain 2 are incorporated as full sequence, and optionally as N-masking

**5)** The SNP information of strain and strain 2 relative to the reference genome build are compared, and a new Ref/SNP annotation is constructed whereby the new Ref/SNP information will be Strain/Strain2 (and no longer the standard reference genome strain Black6 (C57BL/6J)). Later on, genome 1 will correspond to 'strain', genome 2 to 'strain2'.

**6)** The full genome sequence given with `--strain <name>` is read into memory, and the high-confidence SNP positions between Strain and Strain2 are incorporated as full sequence, and optionally as N-masking

The resulting `.fa` files are ready to be indexed with your favourite aligner. Proven and tested aligners include Bowtie2, Tophat, STAR, HISAT2, HiCUP and Bismark. Please note that BWA does not support alignments to N-masked genome (see [here](../SNPsplit/specific_comments.md#bwa)).

Both the SNP filtering and the genome preparation write out report files for record keeping.

## Filtering and processing high confidence SNPs from the VCF file

This section describes in more detail the process of how high confidence SNPs are extracted from the VCF file **mgp_REL2021_snps.vcf.gz**. We are first going to paste the start of the VCF file and explain in more detail later the individual steps taken by the SNPsplit genome preparation for the strain CAST_EiJ as an example; the PDF version of the user guide had relevant information in the header lines or the variant data itself are marked in dark red, but unfotunately this is not supported in this markdown version, so please look carefully... . This should help you adapt the process to other genomes/VCF files should you wish to do so.

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.13+htslib-1.13
##bcftoolsCommand=mpileup -f Mus_musculus.GRCm39.dna.toplevel.fa.gz -b samples -g 10 -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD -E -Q 0 -pm 3 -F 0.25 -d 500
##reference=Mus_musculus.GRCm39.dna.toplevel.fa.gz
##contig=<ID=1,length=195154279>
##contig=<ID=2,length=181755017>
##contig=<ID=3,length=159745316>
##contig=<ID=4,length=156860686>
##contig=<ID=5,length=151758149>
##contig=<ID=6,length=149588044>
##contig=<ID=7,length=144995196>
##contig=<ID=8,length=130127694>
##contig=<ID=9,length=124359700>
##contig=<ID=10,length=130530862>
##contig=<ID=11,length=121973369>
##contig=<ID=12,length=120092757>
##contig=<ID=13,length=120883175>
##contig=<ID=14,length=125139656>
##contig=<ID=15,length=104073951>
##contig=<ID=16,length=98008968>
##contig=<ID=17,length=95294699>
##contig=<ID=18,length=90720763>
##contig=<ID=19,length=61420004>
##contig=<ID=X,length=169476592>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Total allelic depths (high-quality bases)">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled Genotype Quality">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callCommand=call -mAv -f GQ,GP -p 0.99; Date=Wed Aug 11 21:20:03 2021
##bcftools_normCommand=norm --fasta-ref Mus_musculus.GRCm39.dna.toplevel.fa.gz -m +indels; Date=Fri Aug 13 11:11:49 2021
##FORMAT=<ID=FI,Number=1,Type=Integer,Description="High confidence (1) or low confidence (0) based on soft filtering values">
##FILTER=<ID=LowQual,Description="Low quality variants">
##VEP="v104" time="2021-08-30 23:27:00" cache="mus_musculus/104_GRCm39" ensembl-funcgen=104.59ae779 ensembl-variation=104.6154f8b ensembl=104.1af1dce ensembl-io=104.1d3bb6e assembly="GRCm39" dbSNP="150" gencode="GENCODE M27" regbuild="1
.0" sift="sift"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Am
ino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|SIFT|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">
##bcftools_viewVersion=1.13+htslib-1.13
##bcftools_viewCommand=view -i 'FORMAT/FI[*] = 1' mgp_REL2021_snps.vcf.gz; Date=Sat Dec 18 19:08:09 2021
##bcftools_annotateVersion=1.13+htslib-1.13
##bcftools_annotateCommand=annotate -x INFO/VDB,INFO/SGB,INFO/RPBZ,INFO/MQBZ,INFO/MQBZ,INFO/MQSBZ,INFO/BQBZ,INFO/SCBZ,INFO/FS,INFO/MQOF,INFO/AC,INFO/AN,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,FORMAT/GP; Date=Sat Dec 18 19:08:09 2021
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##bcftools_viewCommand=view -a -Oz -o final_mgp_REL2021_snps.vcf.gz; Date=Sat Dec 18 19:08:09 2021
##bcftools_annotateCommand=annotate -x INFO/MQ0F -Oz -o final_mgp_REL2021_snps.vcf.gz mgp_REL2021_snps.vcf.gz; Date=Mon Dec 20 07:12:23 2021
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	129P2_OlaHsd	129S1_SvImJ	129S5SvEvBrd	A_J	AKR_J	B10.RIII	BALB_cByJ	BALB_cJ	BTBR_T+_Itpr3tf_J	BUB_BnJ	C3H_HeH	C3H_HeJ	C57BL_10J	C57BL_10SnJ	C57BL_6NJ	C57BR_cdJ	C57L_J	C58_J	CAST_EiJ	CBA_J	CE_J	CZECHII_EiJ	DBA_1J	DBA_2J	FVB_NJ	I_LnJ	JF1_MsJ	KK_HiJ	LEWES_EiJ	LG_J	LP_J	MAMy_J	MOLF_EiJ	NOD_ShiLtJ	NON_LtJ	NZB_B1NJ	NZO_HlLtJ	NZW_LacJ	PL_J	PWK_PhJ	QSi3	QSi5	RF_J	RIIIS_J	SEA_GnJ	SJL_J	SM_J	SPRET_EiJ	ST_bJ	SWR_J	WSB_EiJ	ZALENDE_EiJ
1	3050050	.	C	G	186.728	PASS	DP=3533;AD=3488,42;DP4=3041,447,44,1;MQ=60;CSQ=G|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||,A|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||;AC=2;AN=104	GT:PL:DP:AD:GQ:FI	0/0:0,57,168:19:19,0:71:1	0/0:0,108,255:36:36,0:122:1	0/0:0,30,154:10:10,0:44:1	0/0:0,199,255:66:66,0:127:1	0/0:0,84,222:28:28,0:98:1	0/0:0,208,255:69:69,0:127:1	0/0:0,163,255:54:54,
0:127:1	0/0:0,223,255:74:74,0:127:1	0/0:0,138,255:46:46,0:127:1	0/0:0,93,254:31:31,0:107:1	0/0:0,9,85:3:3,0:23:0	0/0:0,111,242:37:37,0:125:1	0/0:0,120,255:40:40,0:127:1	0/0:0,154,233:51:51,0:127:1	0/0:0,223,25
5:74:74,0:127:1	0/0:0,255,255:171:171,0:127:0	0/0:0,255,255:213:212,0:127:0	0/0:0,181,255:60:60,0:127:1	0/0:0,163,255:54:54,0:127:1	0/0:0,96,255:32:32,0:110:1	0/0:0,190,255:63:63,0:127:1	0/0:0,255,255:500:499,0:127:
0	0/0:0,93,255:31:31,0:107:1	0/0:0,96,255:32:32,0:110:1	0/0:0,51,233:17:17,0:65:1	0/0:0,90,255:30:30,0:104:1	0/0:0,255,255:220:220,0:127:0	0/0:0,255,255:106:106,0:127:1	0/0:0,196,255:65:65,0:127:0	0/0:
0,18,178:6:6,0:32:1	0/0:0,93,231:31:31,0:107:1	0/0:0,69,255:23:23,0:83:1	0/0:0,255,255:208:208,0:127:0	0/0:0,96,255:32:32,0:110:1	0/0:0,81,255:27:27,0:95:1	0/0:0,111,255:37:37,0:125:1	0/0:0,199,255:66:66,
0:127:1	0/0:0,144,255:49:48,0:127:1	0/0:0,57,255:19:19,0:71:1	0/0:0,255,255:366:366,0:127:0	0/0:0,93,255:31:31,0:107:1	0/0:0,78,255:26:26,0:92:1	0/0:0,69,255:23:23,0:83:1	0/0:0,187,255:62:62,0:127:1	0/0:
0,166,255:55:55,0:127:1	0/0:0,54,224:18:18,0:68:1	0/0:0,12,122:4:4,0:26:0	1/1:252,126,0:42:0,42:105:1	0/0:0,151,255:50:50,0:127:1	0/0:0,114,255:38:38,0:127:1	0/0:0,226,255:75:75,0:127:1	0/0:0,39,201:13:13,0:53:1
1	3050069	.	C	T	5422.31	PASS	DP=3674;AD=2868,796;DP4=2447,421,659,147;MQ=60;CSQ=T|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||,A|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||,G|interge
nic_variant|MODIFIER|||||||||||||||||||SNV||||||||;AC=45;AN=104	GT:PL:DP:AD:GQ:FI	1/1:198,72,0:24:0,24:67:1	1/1:255,120,0:40:0,40:115:1	1/1:149,27,0:9:0,9:22:1	0/0:0,208,255:69:69,0:127:1	1/1:253,90,0:35:1,34:85:1	0/0:0,226,255:75:75,0:127:1	0/0:0,166,255:55:55,0:127:1	0/0:0,199,255:66:66,0:127:1	1/1:255,160,0:53:0,53:127:1	1/1:255,96,0:32:0,32:91:1	0/0:0,12,119:4:4,0:10:0	0/0:0,129,255:43:43,0:127:1	0/0:0,123,255:41:41,
0:121:1	0/0:0,148,255:50:49,0:127:1	0/0:0,223,255:74:74,0:127:1	0/1:200,0,255:175:138,37:127:0	0/1:104,0,255:235:190,44:105:0	0/0:0,163,255:54:54,0:127:1	0/1:179,0,226:55:28,27:127:0	0/0:0,93,255:31:31,0:91:1	0/1:
255,0,255:65:35,30:127:0	0/0:0,255,255:493:467,26:127:0	0/0:0,96,255:32:32,0:94:1	0/0:0,105,255:35:35,0:103:1	1/1:234,72,0:24:0,24:67:1	1/1:255,90,0:30:0,30:85:1	0/0:0,255,255:228:219,8:127:0	0/0:0,255,25
5:109:109,0:127:1	0/0:0,181,255:61:60,0:127:0	0/0:0,24,170:9:8,0:22:1	1/1:239,96,0:32:0,32:91:1	1/1:252,63,0:22:0,21:58:1	0/0:0,255,255:228:226,1:127:0	1/1:255,105,0:35:0,35:100:1	1/1:255,84,0:28:0,28:79:1	0/1:255,0,178:36:17,19:127:0	0/0:0,208,255:70:69,0:127:1	0/0:0,138,255:46:46,0:127:1	1/1:255,75,0:25:0,25:70:1	0/0:0,255,255:387:368,18:127:0	1/1:255,123,0:41:0,41:118:1	1/1:255,87,0:29:0,29:82:1	0/0:0,81,255
:27:27,0:79:1	0/0:0,205,255:68:68,0:127:1	0/0:0,172,255:57:57,0:127:1	1/1:222,57,0:19:0,19:52:1	1/1:142,18,0:6:0,6:13:0	0/0:0,114,255:38:38,0:112:1	1/1:255,148,0:49:0,49:127:1	1/1:255,117,0:39:0,39:112:1	0/0:
0,208,255:70:69,0:127:1	1/1:251,48,0:16:0,16:43:1
…
```

## Detecting strains
To detect the available strains in the VCF file as well as to determine the column number of a desired strain in the file we skim through the header lines until we find a line starting with #CHROM. In terms of the strains in the file the next eight fields are irrelevant: 

`#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT`

The following fields should contain the strains listed in the file, here:

```
129S1_SvImJ	129S5SvEvBrd    A_J	AKR_J	B10.RIII	BALB_cByJ	BALB_cJ	BTBR_T+_Itpr3tf_J	BUB_BnJ	C3H_HeH	C3H_HeJ	C57BL_10J	C57BL_10SnJ	C57BL_6NJ	C57BR_cdJ	C57L_J	C58_J	CAST_EiJ	CBA_J	CE_J	CZECHII_EiJ	DBA_1J	DBA_2J	FVB_NJ	I_LnJ	JF1_MsJ	KK_HiJ	LEWES_EiJ	LG_J	LP_J	MAMy_J	MOLF_EiJ	NOD_ShiLtJ	NON_LtJ	NZB_B1NJ	NZO_HlLtJ	NZW_LacJ	PL_J	PWK_PhJ	QSi3	QSi5	RF_J	RIIIS_J	SEA_GnJ	SJL_J	SM_J	SPRET_EiJ	ST_bJ	SWR_J	WSB_EiJ	ZALENDE_EiJ
```

In our example the data for strain *CAST_EiJ* would be found in column 28, or have an index of 27 (the index = column - 1 because index counting starts at 0 and not at 1). 

If #CHROM is **not present** in the VCF file header, the **automated lookup and subsequent steps will fail**.

## Detecting chromosomes
Chromosomes to be processed are detected from the VCF header lines starting with ##contig=<ID=…,… >, e.g. here:

```
##contig=<ID=1,length=195471971>
##contig=<ID=10,length=130694993>
…
```
The IDs here are the chromosomes for which variant calls are available, but not necessarily all chromosomes in the genome have to be present here. If there are no variant calls for say “chrY” then the subsequent genome preparation step will simply not introduce any changes but use the reference sequence for “chrY”.

If the ##contig=… fields are missing from the VCF file, **subsequent steps are probably going to fail**; as a way around this you might get away with defining the chromosome array manually, e.g. like so:

`my @chroms = (1..19,’X’,’Y’,’MT’);`


## Dealing with Variant Calls
For variant calls the SNPsplit genome preparation extracts the following information from each line:

```
CHROM  [Col 1]
POS    [Col 2]
REF    [Col 4]
ALT    [Col 5]
STRAIN [Col 25 (for CAST_EiJ)]
```

Which in our case is (three lines only):
```
CHROM	POS	REF	ALT	CAST_EiJ
1	3050050	C	G	1/1:0,163,255:54:54,0:127:0
1	3050076	C	T	0/0:0,214,255:71:71,0:127:0	
1	3050069	C	T	1/1:234,72,0:24:0,24:67:1
```

Now the STRAIN field contains a lot of information which is specified in the FORMAT field and further ‘explained’ in the VCF header. The format is: 
```
GT:PL:DP:AD:GQ:FI
```

I am not going into the details about all these FORMAT tags (feel free to browse the header section above), but suffice it to say that the SNPsplit genome preparation only cares about the GT (=GENOTYPE) and FI (=FILTER) entry which are defined as:


```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
```

### GT (Genotype)

The GT string can be one of the following:
```
'.'    = no genotype call was made
'0/0'  = genotype is the same as the reference genome
'1/1'  = homozygous alternative allele; can also be '2/2', '3/3', etc. if more than one alternative allele is present.
'0/1'  = heterozygous genotype; can also be '1/2', '0/2', etc.
```

For the mouse strains we are working with here we only accept homozygous alternative alleles, so '1/1', ‘2/2’ or ‘3/3’. Any other combination is recorded and mentioned in the final report but not included in the SNP file.

### FI (Filter)

```
##FORMAT=<ID=FI,Number=1,Type=Integer,Description="High confidence (1) or low confidence (0) based on soft filtering values">
```

Next we are looking at the FILTER field. If a SNP variant call passed all filters for a given strain it will PASS or get a value of 1. Else it will have a value of 0. In the example above the first position had a genotype of **1/1** but did not pass the Filter for the CAST_EiJ strain (FI value = 0), so this position won’t be included:

```
1	3050050	C	G	1/1:0,163,255:54:54,0:127:0	
```

The second position did not indicate that there was a SNP (**0/0**), so the positions wouldn’t be considered for a single-hybrid genome anyway. However, the filter call for that position was uncertain (FI value = 0), so this position would also not be considered for a SNP position in --dual_hybrid mode even if the other strain had a high confidence SNP at this position (this behaviour was introduced in v0.3.2).

```
1	3050076	C	T	0/0:0,214,255:71:71,0:127:0
```

The third position finally is a high-confidence homozygous SNP which would be included for both single and dual hybrid genomes:

```
1	3050069	C	T	1/1:234,72,0:24:0,24:67:1
```
	
## Clearly defined clean genotype

And finally, variants are required to have a clearly defined genotype, e.g. a made-up position:
```
12	45630185	A	T/C
```

would not be included as a valid SNP position because the alternative allele is not clearly defined (unless the genotype would be 2/2 here, in which case the reference would be ‘A’ and the alternative allele would be ‘C’).



If you are dealing with any VCF file other than the one used as an example here you might need to adapt one or more of the issues addressed here to make it work.
