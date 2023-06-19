#!/usr/bin/env bash
# author :zhang jian
# date:2023.1.14
# version:1.0v
# this is mapping result file deal
# change dir
cd $1
# set run parameter
num_threads=10
# create dir 
mkdir -p  7.bam_file
# sam 2 bam
for i in ./6.mapping_result/*.sam
do
file_name=`basename $i .sam`
# sam 2 bam
samtools view -@ $num_threads -b ./6.mapping_result/$file_name.sam -o ./7.bam_file/$file_name.bam
# bam sorted
samtools sort -@ $num_threads ./7.bam_file/$file_name.bam -o ./7.bam_file/$file_name.sort.bam
# sorted bam index
samtools index -@ $num_threads ./7.bam_file/$file_name.sort.bam
done
