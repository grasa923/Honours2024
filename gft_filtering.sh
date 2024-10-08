##Date:18/08/2024
##Author: Sarah Grant
##Description: Creating new .gtf file containing only exons for BRCA1 to allow exon counts across all samples with feature counts


gtf=/STORAGE/Genomes/Gencode/hg38/gencode.v44.annotation.gtf

## use this to identify ENSE number for exon1
awk '$3 == "exon" && $10 == "\"ENSG00000012048.25\";" && $22 == "1;" {print$0}' $gtf  | grep MANE_Select  |less -S

## make a MANE_Select BRCA1 annotation file
awk '$10 == "\"ENSG00000012048.25\";" {print$0}' $gtf  | grep MANE_Select  > BRCA1_MANE.gtf


featureCounts -T 5 -g 'exon_id' -p --countReadPairs -a BRCA1_MANE.gtf -o ExonCounts.txt align/*.bam