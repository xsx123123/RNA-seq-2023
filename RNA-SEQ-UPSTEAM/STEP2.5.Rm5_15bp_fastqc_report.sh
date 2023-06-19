#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is 2nd fastqc scripts for cutadapter read

# change dir
cd $1 
# create dir
mkdir -p ./2.fastqc_report/rm_5_15bp_fastqc

# use fastqc qc  fastq.gz  file
fastqc -t 10 \
       -q \
       -o ./2.fastqc_report/rm_5_15bp_fastqc \
       ./4.rm_15bp/*.fastq.gz
# meige fastqc report
multiqc ./2.fastqc_report/rm_5_15bp_fastqc \
        -o ./2.fastqc_report/rm_5_15bp_fastqc
# how to use this scripts
# 2.fastqc_report.sh <fastq.gz file path> 