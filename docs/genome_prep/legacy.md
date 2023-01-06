# Legacy (v5/GRCm38)

!!!  warning
    This document references the now outdated mouse genome build GRCm38 and the equally outdated v5 annotation, please consider it as a legacy version only. We do not recommend or support this old version any longer, please go here for the most recent [genome/annotation version](./genome_preparation.md).

`SNPsplit_genome_preparation` is designed to read in a variant call file from the Mouse Genomes Project (e.g. this [v5 file, GRCm38](ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz)) and generate new genome versions where the strain SNPs are either incorporated into the new genome (full sequence) or masked by the ambiguity nucleobase `N` (**N-masking**).

`SNPsplit_genome_preparation` may be run in two different modes:

#### Single strain mode:
**1)** The VCF file is read and filtered for high-confidence SNPs in the strain specified with   strain <name>
**2)** The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated either as N-masking (default), or full sequence (option `--full_sequence`)

#### Dual strain mode:
**1)** The VCF file is read and filtered for high-confidence SNPs in the strain specified with `--strain <name>`
**2)** The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated as full sequence and optionally as N-masking
**3)** The VCF file is read one more time and filtered for high-confidence SNPs in strain 2 specified with `--strain2 <name>`
**4)** The filtered high-confidence SNP positions of strain 2 are incorporated as full sequence, and optionally as N-masking
**5)** The SNP information of strain and strain 2 relative to the reference genome build are compared, and a new Ref/SNP annotation is constructed whereby the new Ref/SNP information will be Strain/Strain2 (and no longer the standard reference genome strain Black6 (C57BL/6J))
**6)** The full genome sequence given with `--strain <name>` is read into memory, and the high-confidence SNP positions between Strain and Strain2 are incorporated as full sequence, and optionally as N-masking

The resulting `.fa` files are ready to be indexed with your favourite aligner. Proved and tested aligners include Bowtie2, Tophat, STAR, HISAT2, HiCUP and Bismark. Please note that STAR and HISAT2 may require you to disable soft-clipping, please see above for more details.

Both the SNP filtering and the genome preparation write out report files for record keeping.

## Filtering and processing high confidence SNPs from the VCF file

This section describes in more detail the process of how high confidence SNPs are extracted from the sample VCF file **mgp.v5.merged.snps_all.dbSNP142.vcf.gz**. We are first going to paste the start of the VCF file and explain in more detail later the individual steps taken by the SNPsplit genome preparation for the strain CAST_EiJ as an example strain; the PDF version of the user guide had relevant information in the header lines or the variant data itself are marked in dark red, but unfotunately this is not supported in this markdown version, so please look carefully... . This should help you adapt the process to other genomes/VCF files should you wish to do so.

```
VCF file: mgp.v5.merged.snps_all.dbSNP142.vcf.gz
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##samtoolsVersion=1.1+htslib-1.1
##bcftools_callVersion=1.1+htslib-1.1
##reference=ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa
##source_20141009.1=vcf-annotate(r953) -f +/D=50/d=5/q=20/w=2/a=5/ (chromosomes=1-19,X,Y: C3H_HeH)
##source_20141009.1=vcf-annotate(r953) -f +/D=100/d=5/q=20/w=2/a=5/ (chromosomes=1-19,X,Y: 129S5SvEvBrd,ZALENDE_EiJ,LEWES_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=200/d=5/q=20/w=2/a=5/ (chromosomes=1-19,X,Y: C57BL_10J)
##source_20141009.1=vcf-annotate(r953) -f +/D=250/d=5/q=20/w=2/a=5/ (chromosomes=1-19,X,Y: 129P2_OlaHsd,A_J,CAST_EiJ,LP_J,PWK_PhJ,WSB_EiJ,BUB_BnJ,DBA_1J,I_LnJ,MOLF_EiJ,NZB_B1NJ,SEA_GnJ,RF_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=300/d=5/q=20/w=2/a=5/ (chromosomes=1-19,X,Y: AKR_J,BALB_cJ,C3H_HeJ,C57BL_6NJ,CBA_J,DBA_2J,C57BR_cdJ,C58_J,NZW_LacJ,C57L_J,KK_HiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=350/d=5/q=20/w=2/a=5/ (chromosomes=1-19,X,Y: 129S1_SvImJ,FVB_NJ,NOD_ShiLtJ,NZO_HlLtJ,SPRET_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=400/d=5/q=20/w=2/a=5/ (chromosomes=1-19,X,Y: BTBR_T+_Itpr3tf_J,ST_bJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=350/d=5/q=20/w=2/a=5/ (chromosome=MT: LEWES_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=650/d=5/q=20/w=2/a=5/ (chromosome=MT: I_LnJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=850/d=5/q=20/w=2/a=5/ (chromosome=MT: BUB_BnJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=1000/d=5/q=20/w=2/a=5/ (chromosome=MT: SEA_GnJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=1200/d=5/q=20/w=2/a=5/ (chromosome=MT: C57BL_10J)
##source_20141009.1=vcf-annotate(r953) -f +/D=1300/d=5/q=20/w=2/a=5/ (chromosome=MT: RF_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=1450/d=5/q=20/w=2/a=5/ (chromosome=MT: KK_HiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=1550/d=5/q=20/w=2/a=5/ (chromosome=MT: NZB_B1NJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=1800/d=5/q=20/w=2/a=5/ (chromosome=MT: ZALENDE_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=2200/d=5/q=20/w=2/a=5/ (chromosome=MT: C3H_HeH)
##source_20141009.1=vcf-annotate(r953) -f +/D=2300/d=5/q=20/w=2/a=5/ (chromosome=MT: ST_bJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=3050/d=5/q=20/w=2/a=5/ (chromosome=MT: C57L_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=3650/d=5/q=20/w=2/a=5/ (chromosome=MT: MOLF_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=4250/d=5/q=20/w=2/a=5/ (chromosome=MT: NZW_LacJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=4650/d=5/q=20/w=2/a=5/ (chromosome=MT: C3H_HeJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=5300/d=5/q=20/w=2/a=5/ (chromosome=MT: C57BR_cdJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=6150/d=5/q=20/w=2/a=5/ (chromosome=MT: DBA_1J)
##source_20141009.1=vcf-annotate(r953) -f +/D=6200/d=5/q=20/w=2/a=5/ (chromosome=MT: 129S5SvEvBrd,C58_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=6650/d=5/q=20/w=2/a=5/ (chromosome=MT: C57BL_6NJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=6700/d=5/q=20/w=2/a=5/ (chromosome=MT: WSB_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=6900/d=5/q=20/w=2/a=5/ (chromosome=MT: CBA_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=7100/d=5/q=20/w=2/a=5/ (chromosome=MT: BALB_cJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=7450/d=5/q=20/w=2/a=5/ (chromosome=MT: NZO_HlLtJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=7650/d=5/q=20/w=2/a=5/ (chromosome=MT: PWK_PhJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=8300/d=5/q=20/w=2/a=5/ (chromosome=MT: SPRET_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=8850/d=5/q=20/w=2/a=5/ (chromosome=MT: CAST_EiJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=9550/d=5/q=20/w=2/a=5/ (chromosome=MT: 129S1_SvImJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=9600/d=5/q=20/w=2/a=5/ (chromosome=MT: LP_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=9850/d=5/q=20/w=2/a=5/ (chromosome=MT: DBA_2J)
##source_20141009.1=vcf-annotate(r953) -f +/D=10200/d=5/q=20/w=2/a=5/ (chromosome=MT: BTBR_T__Itpr3tf_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=11300/d=5/q=20/w=2/a=5/ (chromosome=MT: AKR_J)
##source_20141009.1=vcf-annotate(r953) -f +/D=11550/d=5/q=20/w=2/a=5/ (chromosome=MT: NOD_ShiLtJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=11700/d=5/q=20/w=2/a=5/ (chromosome=MT: 129P2_OlaHsd)
##source_20141009.1=vcf-annotate(r953) -f +/D=11750/d=5/q=20/w=2/a=5/ (chromosome=MT: FVB_NJ)
##source_20141009.1=vcf-annotate(r953) -f +/D=11800/d=5/q=20/w=2/a=5/ (chromosome=MT: A_J)
##contig=<ID=1,length=195471971>
##contig=<ID=10,length=130694993>
##contig=<ID=11,length=122082543>
##contig=<ID=12,length=120129022>
##contig=<ID=13,length=120421639>
##contig=<ID=14,length=124902244>
##contig=<ID=15,length=104043685>
##contig=<ID=16,length=98207768>
##contig=<ID=17,length=94987271>
##contig=<ID=18,length=90702639>
##contig=<ID=19,length=61431566>
##contig=<ID=2,length=182113224>
##contig=<ID=3,length=160039680>
##contig=<ID=4,length=156508116>
##contig=<ID=5,length=151834684>
##contig=<ID=6,length=149736546>
##contig=<ID=7,length=145441459>
##contig=<ID=8,length=129401213>
##contig=<ID=9,length=124595110>
##contig=<ID=X,length=171031299>
##contig=<ID=Y,length=91744698>
##contig=<ID=MT,length=16299>
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##QUAL=<ID=QUAL,Number=1,Type=Float,Description="The highest QUAL value for a variant location from any of the samples">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Total Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type from Ensembl 75 as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled genotype posterior probabilities">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of high-quality non-reference bases">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##FORMAT=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##FORMAT=<ID=FI,Number=1,Type=Integer,Description="Whether a sample was a Pass(1) or fail (0) based on FILTER values">
##FILTER=<ID=StrandBias,Description="Min P-value for strand bias (INFO/PV4) [0.0001]">
##FILTER=<ID=EndDistBias,Description="Min P-value for end distance bias (INFO/PV4) [0.0001]">
##FILTER=<ID=MaxDP,Description="Maximum read depth (INFO/DP or INFO/DP4) []">
##FILTER=<ID=BaseQualBias,Description="Min P-value for baseQ bias (INFO/PV4) [0]">
##FILTER=<ID=MinMQ,Description="Minimum RMS mapping quality for SNPs (INFO/MQ) [20]">
##FILTER=<ID=MinAB,Description="Minimum number of alternate bases (INFO/DP4) [5]">
##FILTER=<ID=Qual,Description="Minimum value of the QUAL field [10]">
##FILTER=<ID=VDB,Description="Minimum Variant Distance Bias (INFO/VDB) [0]">
##FILTER=<ID=GapWin,Description="Window size for filtering adjacent gaps [3]">
##FILTER=<ID=MapQualBias,Description="Min P-value for mapQ bias (INFO/PV4) [0]">
##FILTER=<ID=SnpGap,Description="SNP within INT bp around a gap to be filtered [2]">
##FILTER=<ID=RefN,Description="Reference base is N []">
##FILTER=<ID=MinDP,Description="Minimum read depth (INFO/DP or INFO/DP4) [5]">
##FILTER=<ID=Het,Description="Genotype call is heterozygous (low quality) []">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  129P2_OlaHsd    129S1_SvImJ     129S5SvEvBrd    AKR_J   A_J     BALB_cJ BTBR_T+_Itpr3tf_J       BUB_BnJ C3H_HeH C3H_HeJ C57BL_10J       C57BL_6NJ       C57BR_cdJ   C57L_J   C58_J   CAST_EiJ        CBA_J   DBA_1J  DBA_2J  FVB_NJ  I_LnJ   KK_HiJ  LEWES_EiJ       LP_J    MOLF_EiJ        NOD_ShiLtJ      NZB_B1NJ        NZO_HlLtJ       NZW_LacJ        PWK_PhJ RF_J    SEA_GnJ SPRET_EiJ       ST_bJ   WSB_EiJ      ZALENDE_EiJ
1       3000023 .       C       A       153     MinDP;MinAB;Qual;Het    DP=170;DP4=2,0,168,0;CSQ=A||||intergenic_variant||||||||        GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  1/1:20:8:0:133,20,0:109,11,0:2:29:7:1,0,7,0:0:-0.6364
26:.:1  1/1:22:6:0.166667:152,22,0:137,18,0:2:36:6:0,0,6,0:0:-0.616816:.:1      1/1:11:4:0:70,11,0:51,4,0:2:24:3:1,0,3,0:0:-0.511536:.:0        1/1:15:4:0.25:80,15,0:68,12,0:2:25:4:0,0,4,0:0:-0.556411:.:0    1/1:11:3:0:72,11,0:63,9,0:2:3
4:3:0,0,3,0:0:-0.511536:.:0     0/1:3:1:0:40,3,3:37,3,0:2:40:1:0,0,1,0:0:-0.379885:.:0  1/1:36:10:0:194,36,0:174,30,0:2:38:10:0,0,10,0:0:-0.670168:.:1  1/1:22:6:0.333333:111,22,0:96,18,0:2:20:6:0,0,6,0:0:-0.616816:.:1       0/1:3:1:0:35,
3,3:32,3,0:2:40:1:0,0,1,0:0:-0.379885:.:0       1/1:19:5:0:135,19,0:121,15,0:2:34:5:0,0,5,0:0:-0.590765:.:1     1/1:15:4:0:106,15,0:94,12,0:2:29:4:0,0,4,0:0:-0.556411:.:0      1/1:11:3:0:83,11,0:74,9,0:2:28:3:0,0,3,0:0:-0.511536:.:0    1/1:33:9:0:199,33,0:180,27,0:2:42:9:0,0,9,0:0:-0.662043:.:1      1/1:6:2:0:56,6,0:50,6,0:2:27:2:0,0,2,0:0:-0.453602:.:0  1/1:15:4:0.5:75,15,0:63,12,0:2:25:4:0,0,4,0:0:-0.556411:.:0     1/1:15:4:0:79,15,0:67,12,0:2:24:4:0,0,4,0:0:-0.556411
:.:0    1/1:22:6:0:133,22,0:118,18,0:2:31:6:0,0,6,0:0:-0.616816:.:1     1/1:15:4:0:108,15,0:96,12,0:2:35:4:0,0,4,0:0:-0.556411:.:0      1/1:15:4:0.25:88,15,0:76,12,0:2:31:4:0,0,4,0:0:-0.556411:.:0    0/0:.:1:0:.,.,.:.,.,.:2:40:1:0,0,1,0:
0:-0.379885:.:0 1/1:11:3:0:79,11,0:70,9,0:2:31:3:0,0,3,0:0:-0.511536:.:0        1/1:6:2:0.5:48,6,0:42,6,0:2:20:2:0,0,2,0:0:-0.453602:.:0        1/1:15:4:0.5:70,15,0:58,12,0:2:22:4:0,0,4,0:0:-0.556411:.:0     1/1:6:2:0:69,6,0:63,6,0:2:40:
2:0,0,2,0:0:-0.453602:.:0       1/1:19:5:0.2:94,19,0:80,15,0:2:24:5:0,0,5,0:0:-0.590765:.:1     1/1:15:4:0:87,15,0:75,12,0:2:30:4:0,0,4,0:0:-0.556411:.:0       1/1:26:7:0.142857:150,26,0:134,21,0:2:31:7:0,0,7,0:0:-0.636426:.:1      1/1:1
9:5:0:132,19,0:118,15,0:2:48:5:0,0,5,0:0:-0.590765:.:1  1/1:30:8:0:171,30,0:153,24,0:2:37:8:0,0,8,0:0:-0.651104:.:1     1/1:11:3:0.333333:72,11,0:63,9,0:2:33:3:0,0,3,0:0:-0.511536:.:0 1/1:19:5:0:137,19,0:123,15,0:2:38:5:0,0,5,0:0:-0.5907
65:.:1  1/1:22:6:0.166667:137,22,0:122,18,0:2:42:6:0,0,6,0:0:-0.616816:.:1      1/1:33:9:0:140,33,0:121,27,0:2:23:9:0,0,9,0:0:-0.662043:.:1     1/1:33:9:0.111111:170,33,0:151,27,0:2:29:9:0,0,9,0:0:-0.662043:.:1      1/1:26:7:0:140,26,0:1
24,21,0:2:31:7:0,0,7,0:0:-0.636426:.:1  1/1:6:2:0:48,6,0:42,6,0:2:21:2:0,0,2,0:0:-0.453602:.:0
1       3000126 rs580370473     G       T       184     MinDP;MinMQ;Het;MinAB;Qual      DP=417;DP4=93,1,210,113;CSQ=T||||intergenic_variant||||||||     GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI  1/1:4:8:0:89,4,1:75,0,1:2:36:5:3,0,4,
1:0:-0.590765:.:1       1/1:21:12:0:152,21,0:128,12,0:2:55:9:3,0,7,2:0:-0.662043:.:1    0/0:.:2:0:.,.,.:.,.,.:2:16:0:2,0,0,0:0:.:.:0    0/1:9:11:0:91,0,9:82,0,11:2:48:9:2,0,6,3:0:-0.662043:.:0        1/1:6:5:0:67,6,1:51,1,0:2:44:3:2,0,1,
2:4:-0.511536:.:0       1/1:5:10:0:97,5,1:82,0,0:2:56:8:2,0,7,1:0:-0.651104:.:1 1/1:31:28:0:173,31,0:151,23,0:2:50:21:7,0,6,15:28:-0.692352:.:1 1/1:18:8:0:135,18,0:110,9,0:2:56:6:2,0,3,3:3:-0.616816:.:1      1/1:11:3:0:55,11,0:46,9,0:2:3
7:3:0,0,1,2:0:-0.511536:.:0     1/1:50:14:0:111,50,0:89,42,0:2:41:14:0,0,14,0:0:-0.686358:.:1   1/1:23:7:0:115,23,0:86,12,0:2:38:6:1,0,1,5:0:-0.616816:.:1      1/1:6:8:0:66,6,1:52,2,0:2:49:5:3,0,5,0:0:-0.590765:.:1  1/1:60:17:0.0588235:2
17,60,0:193,51,0:2:52:17:0,0,12,5:0:-0.690438:.:1       1/1:75:33:0:244,75,0:211,62,0:2:52:30:3,0,17,13:6:-0.693097:.:1 1/1:47:13:0:217,47,0:195,39,0:2:49:13:0,0,4,9:0:-0.683931:.:1   0/0:.:5:0:.,.,.:.,.,.:2:47:4:1,0,4,0:0:-0.556411:.:01/1:43:12:0:117,43,0:96,36,0:2:47:12:0,0,12,0:0:-0.680642:.:1    1/1:29:8:0:189,29,0:155,15,0:2:53:7:1,0,3,4:0:-0.636426:.:1     1/1:15:4:0:73,15,0:61,12,0:2:59:4:0,0,4,0:0:-0.556411:.:0       0/1:5:6:0:76,1,5:66,0,7:2:50:4:2,0,4,0:0:-0.5
56411:.:0       1/1:13:9:0:97,13,0:78,7,0:2:45:6:3,0,4,2:3:-0.616816:.:1        1/1:37:15:0.0666667:195,37,0:165,25,0:2:52:13:2,0,7,6:3:-0.683931:.:1   1/1:11:21:0.666667:38,11,0:35,12,0:2:9:6:15,0,3,3:18:-0.616816:.:0      1/1:19:6:0:92
,19,0:67,10,0:2:47:5:1,0,5,0:0:-0.590765:.:1    1/1:10:30:0.4:101,10,0:88,7,0:2:24:14:16,0,0,14:82:-0.686358:.:1        0/1:21:7:0:51,0,21:46,0,21:2:48:5:2,0,4,1:0:-0.590765:.:0       1/1:14:7:0.142857:80,14,0:60,7,0:2:44:5:1,1,4,1:0:-0.
590765:.:1      1/1:5:14:0:101,5,1:86,0,0:2:52:9:5,0,9,0:0:-0.662043:.:1        1/1:33:9:0:196,33,0:177,27,0:2:53:9:0,0,6,3:0:-0.662043:.:1     1/1:26:7:0:91,26,0:75,21,0:2:50:7:0,0,7,0:0:-0.636426:.:1       1/1:44:24:0:138,44,0:119,38,0
:2:45:23:1,0,19,4:0:-0.692717:.:1       1/1:17:8:0:101,17,0:79,9,0:2:53:6:2,0,4,2:0:-0.616816:.:1       0/1:20:11:0.0909091:72,0,20:65,0,21:2:45:7:4,0,7,0:0:-0.636426:.:0      1/1:33:23:0:183,33,0:159,24,0:2:50:18:5,0,6,12:19:-0.691153:.
:1      1/1:11:7:0:66,11,0:49,6,0:2:54:5:2,0,5,0:0:-0.590765:.:1        1/1:19:5:0:62,19,0:48,15,0:2:54:5:0,0,5,0:0:-0.590765:.:1
…
…
```

### Detecting strains
To detect the available strains in the VCF file as well as to determine the column number of a desired strain in the file we skim through the header lines until we find a line starting with #CHROM. In terms of the strains in the file the next eight fields are irrelevant: 

`#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT`

The following fields should contain the strains listed in the file, here: 
`129P2_OlaHsd	129S1_SvImJ	129S5SvEvBrd	AKR_J	A_J	BALB_cJ	BTBR_T+_Itpr3tf_J	BUB_BnJ	C3H_HeH	C3H_HeJ	C57BL_10J	C57BL_6NJ	C57BR_cdJ	C57L_J	C58_J	CAST_EiJ	CBA_J	DBA_1J …`

In our example the data for strain *CAST_EiJ* would be found in column 25, or have an index of 24 (the index = column - 1 because index counting starts at 0 and not at 1). 

If #CHROM is **not present** in the VCF file header, the **automated lookup and subsequent steps will fail**.

### Detecting chromosomes
Chromosomes to be processed are detected from the VCF header lines starting with ##contig=<ID=…,… >, e.g. here:

```
##contig=<ID=1,length=195471971>
##contig=<ID=10,length=130694993>
…
```
The IDs here are the chromosomes for which variant calls are available, but not necessarily all chromosomes in the genome have to be present here. If there are no variant calls for say “chrY” then the subsequent genome preparation step will simply not introduce any changes but use the reference sequence for “chrY”.

If the ##contig=… fields are missing from the VCF file, **subsequent steps are probably going to fail**; as a way around this you might get away with defining the chromosome array manually, e.g. like so:

`my @chroms = (1..19,’X’,’Y’,’MT’);`


### Dealing with Variant Calls
For variant calls the SNPsplit genome preparation extracts the following information from each line:

```
CHROM  [Col 1]
POS    [Col 2]
REF    [Col 4]
ALT    [Col 5]
STRAIN [Col 25 (for CAST_EiJ)]
```

Which in our case is (first three lines only):
```
CHROM	POS	REF	ALT	CAST_EiJ
1	3000023	C	A	1/1:15:4:0:79,15,0:67,12,0:2:24:4:0,0,4,0:0:-0.556411:.:0
1	3000126	G	T	0/0:.:5:0:.,.,.:.,.,.:2:47:4:1,0,4,0:0:-0.556411:.:0
1	3000185	G	T	1/1:43:12:0:276,43,0:255,36,0:2:54:12:0,0,10,2:0:-0.680642:.:1
```

Now the STRAIN field contains a lot of information which is specified in the FORMAT field and further ‘explained’ in the VCF header. The format is: 
```
GT:GQ:DP:MQ0F:GP:PL:AN:MQ:DV:DP4:SP:SGB:PV4:FI
```

I am not going into the details about all these FORMAT tags (feel free to browse the header section above), but suffice it to say that the SNPsplit genome preparation only cares about the GT (=GENOTYPE) and FI (=FILTER) entry which are defined as:


```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
```

The GT string can be one of the following:
```
'.'    = no genotype call was made
'0/0'  = genotype is the same as the reference genome
'1/1'  = homozygous alternative allele; can also be '2/2', '3/3', etc. if more than one alternative allele is present.
'0/1'  = heterozygous genotype; can also be '1/2', '0/2', etc.
```

For the mouse strains we are working with here we only accept homozygous alternative alleles, so '1/1', ‘2/2’ or ‘3/3’. Any other combination is recorded and mentioned in the final report but not included in the SNP file.

```
##FORMAT=<ID=FI,Number=1,Type=Integer,Description="Whether a sample was a Pass(1) or fail (0) based on FILTER values">
```

Next we are looking at the FILTER field. If a SNP variant call passed all filters for a given strain it will PASS or get a value of 1. Else it will have a value of 0. In the example above the first position had a genotype of **1/1** but did not pass the Filter for the CAST_EiJ strain (FI value = 0), so this position won’t be included:

```
1	3000023	C	A	1/1:15:4:0:79,15,0:67,12,0:2:24:4:0,0,4,0:0:-0.556411:.:0
```

The second position did not indicate that there was a SNP (**0/0**), so the positions wouldn’t be considered for a single-hybrid genome anyway. However, the filter call for that position was uncertain (FI value = 0), so this position would also not be considered for a SNP position in --dual_hybrid mode even if the other strain had a high confidence SNP at this position (this behaviour was introduced in v0.3.2).

```
1	3000126	G	T	0/0:.:5:0:.,.,.:.,.,.:2:47:4:1,0,4,0:0:-0.556411:.:0
```

The third position finally is a high-confidence homozygous SNP which would be included for both single and dual hybrid genomes:

```
1	3000185	G	T	1/1:43:12:0:276,43,0:255,36,0:2:54:12:0,0,10,2:0:-0.680642:.:1
```
	
### Clearly defined clean genotype	
And finally, variants are required to have a clearly defined genotype, e.g. a made-up position:
```
12	45630185	A	T/C
```

would not be included as a valid SNP position because the alternative allele is not clearly defined (unless the genotype would be 2/2 here, in which case the reference would be ‘A’ and the alternative allele would be ‘C’).



If you are dealing with any other VCF file than the one used as an example here you might need to adapt one or more of the issues addressed here to make it work.
