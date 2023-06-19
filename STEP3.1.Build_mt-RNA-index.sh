#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is 2nd fastqc scripts for cutadapter read
# set run parameter
rrna_dir=/data/genome/new_mus/ensembl_mm39/rRNA/mus_rRNA.fa
rrna_index=/data/genome/new_mus/ensembl_mm39/rRNA/mus_rRNA_hisat2_index=
# change dir
cd $rrna_dir
# building index by hisat2
hisat2-build -p 10 mus_rRNA.fa ./mus_rRNA_hisat2_index/index