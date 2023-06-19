#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is build index scripts by hisat2
# create dir
cd  $1
# input and output file dir
genome_fasta=/data/genome/new_mus/ensembl_mm39/genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.chr1-mt.fa
hisat2_index=/data/genome/new_mus/ensembl_mm39/index/hisat2_chr1-mt_index/index

# create star index dir
mkdir $hisat2_index
# create STAR index
hisat2-build -p 10 $genome_fasta $hisat2_index/genome_index
echo "hisat2 build index finish"