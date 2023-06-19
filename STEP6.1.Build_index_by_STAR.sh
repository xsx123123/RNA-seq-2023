#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is build index scripts by STAR
# create dir
# cd  $1
# input and output file dir
genome_fasta=/data/genome/new_mus/ensembl_mm39/genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.chr1-mt.fa
gtf=/data/genome/new_mus/ensembl_mm39/annotation/Mus_musculus.GRCm39.109.gtf
STAR_index=/data/genome/new_mus/ensembl_mm39/index/STAR_mm39_index/mm39_index
num_threads=20
# create star index dir
mkdir $STAR_index
# create STAR index
STAR --genomeSAindexNbases 20 \
     --runThreadN $num_threads \
     --runMode genomeGenerate \
     --genomeDir $STAR_index \
     --genomeFastaFiles $genome_fasta \
     --sjdbGTFfile $gtf

echo "STAR build index finish"