#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is mapping scripts by STAR
# create dir
cd $1
# create dir 
mkdir 6.mapping_result
# set run parameter
num_threads=15
index=/data/genome/new_mus/ensembl_mm39/index/STAR_mm39_index/mm39_index
# mapping
for each in ./5.rm_rRNA/not_mapping_read/*.unalign.fastq.1.gz
do
file_name=`basename $each .unalign.fastq.1.gz`
STAR --outFilterType BySJout \
     --runThreadN $num_threads \
     --outFilterMismatchNmax 2 \
     --genomeDir $index \
     --readFilesIn ./5.rm_rRNA/not_mapping_read/$file_name.unalign.fastq.1.gz ./5.rm_rRNA/not_mapping_read/$file_name.unalign.fastq.2.gz \
     --outFileNamePrefix ./6.mapping_result/$file_name \
     --limitBAMsortRAM 100000000000 \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts \
     --outFilterMultimapNmax 1 \
     --outFilterMatchNmin 20 \
     --alignEndsType Local
done
# alignEndsType Local
# outFilterMismatchNmax 2
# outFilterMultimapNmax 1 or less than 3 is OK