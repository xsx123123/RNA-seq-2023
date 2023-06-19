#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is 1st fastqc scripts for fastq.gz file

# change dir
cd $1 
# create dir
mkdir -p  2.fastqc_report/1st_fastqc_report
# use fastqc qc  fastq.gz  file
fastqc -t 10 \
       -q \
       -o ./2.fastqc_report/1st_fastqc_report \
       ./1.raw_data/*fastq.gz
# meige fastqc report
multiqc ./2.fastqc_report/1st_fastqc_report \
        -o ./2.fastqc_report/1st_fastqc_report
# how to use this scripts
# 2.fastqc_report.sh <fastq.gz file path> 