#!/bin/bash

## General paths
path_genomes=/STORAGE/Genomes/Gencode/hg38
genome=/WORKSPACE/George/genomes/STAR_genomes

## test to see if genome has been indexed - if not make a new one based on gencode v44
if [ -f $genome"/sjdbInfo.txt" ]
then
	echo "$genome director will be used"
else
	echo  STAR --runMode genomeGenerate --genomeDir ./indexes/ \
            --genomeFastaFiles $path_genomes"/GRCh38.primary_assembly.genome.fa.gz" \
            --sjdbGTFfile $path_genomes"/gencode.v44.annotation.gtf" \
            --sjdbOverhang 50 --outFileNamePrefix hg38_star/
	genome=hg38_star/
fi


## Align reads and count them using STAR
for read1 in /HANGER/sarahGrant/RNA-seq/Raw/*/*1.fq.gz
do
        read2=${read1%%1.fq.gz}"2.fq.gz"
	out=${read1##*/}
	out=${out%_1.fq.gz}
	STAR --genomeDir $genome \
		--readFilesIn $read1 $read2 \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode GeneCounts \
		--runThreadN 20 \
		--outFileNamePrefix align/$out

        echo ${out} "Complete"
done

