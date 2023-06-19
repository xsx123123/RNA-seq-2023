#!/usr/bin/env bash
# author :zhang jian
# date:2023.1.14
# version:1.0v
# this is mapping to genome scripts (suite deal dutp data)
# change dir
cd $1
# set run parameter
num_threads=15
index=/data/genome/new_mus/ensembl_mm39/index/hisat2_chr1-mt_index/index
# create dir 
mkdir -p 5.mapping_result/mapping_summary
mkdir -p 5.mapping_result/not_mapping_read
mkdir 6.bam_file
# mapping & convert sam 2 sam
for each in ./4.rm_rRNA/not_mapping_read/*.unalign.fastq.1.gz
do
file_name=`\basename $each .unalign.fastq.1.gz`
hisat2 -t \
       -p $threads \
       -k 1 \
       --rna-strandness RF \
       -x $index \
       -1 ./4.rm_rRNA/not_mapping_read/$file_name.unalign.fastq.1.gz \
       -2 ./4.rm_rRNA/not_mapping_read/$file_name.unalign.fastq.2.gz \
       -S ./5.mapping_result/$file_name.sam \
       --un-conc-gz ./5.mapping_result/not_mapping_read/$file_name.unalign.fastq.gz \
       --summary-file  ./5.mapping_result/mapping_summary/${file_name}_report.txt
echo $file_name "mapping Done!!!"
# use samtools convert sam file
samtools view -@ $num_threads  -b ./5.mapping_result/$file_name.sam -o ./6.bam_file/$file_name.bam
# use samtools sort bam
samtools sort -@ $num_threads ./6.bam_file/$file_name.bam -o ./6.bam_file/$file_name.sort.bam
# build index for bam file
samtools index -@ $num_threads ./6.bam_file/$file_name.sort.bam
done

# -k 1 is reserve every read only one uniquely mapping 
# --rna-strandness RF is dutp building lirbry 