#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is remove rrna scripts by hisat2
# set run parameter
rrna_index=/data/genome/new_mus/ensembl_mm39/rRNA/mus_rRNA_hisat2_index/index
# change dir
cd $1
# create dir
mkdir -p 5.rm_rRNA/mapping_rrna
mkdir -p 5.rm_rRNA/not_mapping_read
mkdir -p 5.rm_rRNA/report
# remove rrna
for each in ./4.rm_15bp/*.R1.fastq.gz
do
file_name=`\basename $each .R1.fastq.gz`
echo "mapping" $file_name "is doing"
hisat2 -t \
       -p 20 \
       --rna-strandness RF \
       -x  $rrna_index \
       -1 ./4.rm_15bp/$file_name.R1.fastq.gz \
       -2 ./4.rm_15bp/$file_name.R2.fastq.gz \
       -S ./5.rm_rRNA/mapping_rrna/$file_name.sam \
       --un-conc-gz ./5.rm_rRNA/not_mapping_read/$file_name.unalign.fastq.gz  \
       --summary-file  ./5.rm_rRNA/report/${file_name}_report.txt
echo $file_name "Done!!!"
done

# convert file 1 name
rename 's/.fastq.1./.1.fastq./' ./5.rm_rRNA/not_mapping_read/*.gz
# convert file 2 name
rename 's/.fastq.2./.2.fastq./' ./5.rm_rRNA/not_mapping_read/*.gz