#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is mapping result QC & build index
# set run parameter
num_threads=20
# change dir
cd $1
# create dir
mkdir -p 6.mapping_result/genome_mapping
mkdir -p 2.fastqc_report/genome_mapping_result
mkdir -p 6.mapping_result/transc_mapping
mkdir -p 2.fastqc_report/transc_mapping_result
# move file
mv 6.mapping_result/*Aligned.sortedByCoord.out.bam 6.mapping_result/genome_mapping
mv 6.mapping_result/*Aligned.toTranscriptome.out.bam 6.mapping_result/transc_mapping
# mapping result summary
fastqc -t $num_threads \
       -q 6.mapping_result/genome_mapping/*.bam \
       -o 2.fastqc_report/genome_mapping_result

fastqc -t $num_threads \
       -q 6.mapping_result/transc_mapping/*.bam \
       -o 2.fastqc_report/transc_mapping_result

# merge fastqc report
multiqc ./2.fastqc_report/genome_mapping_result -o ./2.fastqc_report/genome_mapping_result
multiqc ./2.fastqc_report/transc_mapping_result -o ./2.fastqc_report/transc_mapping_result
# build index
# build index for genome
for i in 6.mapping_result/genome_mapping/*.bam
do
samtools index -@ $num_threads \
         $i
done
# do

