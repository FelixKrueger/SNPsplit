# Examples

## Paired-end report (2x50bp):

```
Input file:					'FVBNJ_Cast.bam'
Writing allele-flagged output file to:	'FVBNJ_Cast.allele_flagged.bam'


Allele-tagging report
=====================
Processed 194564995 read alignments in total
149380724 reads were unassignable (76.78%)
35143075 reads were specific for genome 1 (18.06%)
9860248 reads were specific for genome 2 (5.07%)
118662 reads did not contain one of the expected bases at known SNP positions (0.06%)
180948 contained conflicting allele-specific SNPs (0.09%)


SNP coverage report
===================
N-containing reads:	45276050
non-N:			149262062
total:			194564995
Reads had a deletion of the N-masked position (and were thus dropped):	26883 (0.01%)
Of which had multiple deletions of N-masked positions within the same read:	30

Of valid N containing reads,
N was present in the list of known SNPs:		61087551 (99.99%)
N was not present in the list of SNPs:		4773 (0.01%)

Input file:						'FVBNJ_Cast.allele_flagged.bam'
Writing unassigned reads to:				'FVBNJ_Cast.unassigned.bam'
Writing genome 1-specific reads to:			'FVBNJ_Cast.genome1.bam'
Writing genome 2-specific reads to:			'FVBNJ_Cast.genome2.bam'


Allele-specific paired-end sorting report
=========================================
Read pairs/singletons processed in total:		98215744
	thereof were read pairs:			96349251
	thereof were singletons:			1866493
Reads were unassignable (not overlapping SNPs):	61174812 (62.29%)
	thereof were read pairs:	59662537
	thereof were singletons:	1512275
Reads were specific for genome 1:			28657857 (29.18%)
	thereof were read pairs:	28446094
	thereof were singletons:	211763
Reads were specific for genome 2:			8122687 (8.27%)
	thereof were read pairs:	7985424
	thereof were singletons:	137263
Reads contained conflicting SNP information:		260388 (0.27%)
	thereof were read pairs:	255196
	thereof were singletons:	5192
```
 
## Hi-C report (2x100bp):

```
Input file:					Black6_129S1.bam
Writing allele-flagged output file to:	Black6_129S1.allele_flagged.bam


Allele-tagging report
=====================
Processed 94887256 read alignments in total
59662038 reads were unassignable (62.88%)
19851697 reads were specific for genome 1 (20.92%)
15047281 reads were specific for genome 2 (15.86%)
47261 reads did not contain one of the expected bases at known SNP positions (0.05%)
326240 contained conflicting allele-specific SNPs (0.34%)


SNP coverage report
===================
N-containing reads:	35231977
non-N:			59614777
total:			94887256
Reads had a deletion of the N-masked position (and were thus dropped):	40502 (0.04%)
Of which had multiple deletions of N-masked positions within the same read:	59

Of valid N containing reads,
N was present in the list of known SNPs:	57101748 (99.99%)
N was not present in the list of SNPs:	4211 (0.01%)

Input file:						Black6_129S1.allele_flagged.bam'
Writing unassigned reads to:				Black6_129S1.UA_UA.bam'
Writing genome 1-specific reads to:			Black6_129S1.G1_G1.bam'
Writing genome 2-specific reads to:			Black6_129S1.G2_G2.bam'
Writing G1/UA reads to:				Black6_129S1.G1_UA.bam'
Writing G2/UA reads to:				Black6_129S1.G2_UA.bam'
Writing G1/G2 reads to:				Black6_129S1.G1_G2.bam'


Allele-specific paired-end sorting report
=========================================
Read pairs processed in total:			47443628
Read pairs were unassignable (UA/UA):			18862725 (39.76%)
Read pairs were specific for genome 1 (G1/G1):	3533932 (7.45%)
Read pairs were specific for genome 2 (G2/G2):	2592040 (5.46%)
Read pairs were a mix of G1 and UA:			12306421 (25.94%). Of these,
			were G1/UA: 6018598
			were UA/G1: 6287823
Read pairs were a mix of G2 and UA:			9430675 (19.88%). Of these,
			were G2/UA: 4603429
			were UA/G2: 4827246
Read pairs were a mix of G1 and G2:			395296 (0.83%). Of these,
			were G1/G2: 198330
			were G2/G1: 196966
Read pairs contained conflicting SNP information:	322539 (0.68%)
```


 
## BS-Seq report (2x100bp):

```
Input file:					'129_Cast_bismark_bt2_pe.bam'
Writing allele-flagged output file to:	'129_Cast_bismark_bt2_pe.allele_flagged.bam'


Allele-tagging report
=====================
Processed 162441396 read alignments in total
Reads were unaligned and hence skipped: 0 (0.00%)
109109113 reads were unassignable (67.17%)
30267901 reads were specific for genome 1 (18.63%)
22697499 reads were specific for genome 2 (13.97%)
15807753 reads did not contain one of the expected bases at known SNP positions (9.73%)
366883 contained conflicting allele-specific SNPs (0.23%)


SNP coverage report
===================
SNP annotation file:	../all_Cast_SNPs_129S1_reference.mgp.v4.txt.gz
N-containing reads:	68984287
non-N:			93301360
total:			162441396
Reads had a deletion of the N-masked position (and were thus dropped):	155749 (0.10%)
Of which had multiple deletions of N-masked positions within the same read:	65

Of valid N containing reads,
N was present in the list of known SNPs:	119119643 (99.99%)
Positions were skipped since they involved C>T SNPs:	38464451
N was not present in the list of SNPs:		7517 (0.01%)

Input file:						129_Cast_bismark_bt2_pe.allele_flagged.bam'
Writing unassigned reads to:				129_Cast_bismark_bt2_pe.unassigned.bam'
Writing genome 1-specific reads to:			129_Cast_bismark_bt2_pe.genome1.bam'
Writing genome 2-specific reads to:			129_Cast_bismark_bt2_pe.genome2.bam'


Allele-specific paired-end sorting report
=========================================
Read pairs/singletons processed in total:		81220698
	thereof were read pairs:			81220698
	thereof were singletons:			0
Reads were unassignable (not overlapping SNPs):		40420625 (49.77%)
	thereof were read pairs:	40420625
	thereof were singletons:	0
Reads were specific for genome 1:			23037433 (28.36%)
	thereof were read pairs:	23037433
	thereof were singletons:	0
Reads were specific for genome 2:			17303663 (21.30%)
	thereof were read pairs:	17303663
	thereof were singletons:	0
Reads contained conflicting SNP information:		458977 (0.57%)
	thereof were read pairs:	458977
	thereof were singletons:	0
```