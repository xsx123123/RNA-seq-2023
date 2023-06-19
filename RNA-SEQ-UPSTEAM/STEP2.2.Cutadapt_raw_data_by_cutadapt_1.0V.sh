#!/usr/bin/env bash
# author:zhang jian
# date:2023.1.9
# version:1.0v
# e-maiil:zhangjian199567@outlook.com
# this is cutadapt scripts for fastq.gz file

# create dir
cd  $1
# mkdir cutadapt fastq.gz file
mkdir $1/3.rm_adapter

# cutadapt
for each in ./1.raw_data/*.R1.fastq.gz
do
file_name=`basename $each .R1.fastq.gz`
cutadapt -j 20 \
         -q 30 \
         -m 50 \
         -a AGATCGGAAGAG \
         -A AGATCGGAAGAG \
         -o ./3.rm_adapter/$file_name.R1.fastq.gz \
         -p ./3.rm_adapter/$file_name.R2.fastq.gz \
         -O 5 \
         -e 0.1 \
         ./1.raw_data/$file_name.1.fastq.gz \
         ./1.raw_data/$file_name.2.fastq.gz
done
# 单端测序参数
# -g: #剪切reads 5'端adapter(双端测序第一条read)，加$表示adapter锚定在reads 5'端
# -a: #剪切reads 3'端adapter(双端测序第一条read)，加$表示adapter锚定在reads3'端
# -O MINLENGTH, --overlap=MINLENGTH #adapter与reads最小overlap,才算成功识别; Default: 3
# -m LENGTH, --minimum-length=LENGTH: 根据最短长度筛选reads;Default: 0
# --discard-untrimmed, --trimmed-only #丢掉不包含adapter的reads
# -e ERROR_RATE, --error-rate=ERROR_RATE  #adapter匹配允许的最大错配率（错配/匹配片段长度)；Default: 0.1
# --no-trim: 不修剪adapter，直接输出满足跳进啊的reads
# -u LENGTH, --cut=LENGTH:  #修剪reads 5'/3'端碱基,正数：从开始除移除碱基；负数：从末尾处移除碱基；
# -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff=[5'CUTOFF,]3'CUTOFF: #修剪低质量碱基
# -l LENGTH, --length=LENGTH: #将reads修剪的最终长度
# --trim-n: #修剪reads末端的'N'
# -o FILE, --output=FILE: #输出文件
# --info-file=FILE：每条reads和与reads匹配的adapter的信息
# --too-short-output=FILE: #为reads长度最小值设定阈值筛选reads后，要丢弃的部分输出到文件；长度依据m值设定；   
# --too-long-output=FILE：#为reads长度最大值设定阈值筛选reads后，要丢弃的部分输出到文件；长度依据M值设定； 
# --untrimmed-output=FILE: #将没有adapter未做修剪的reads输出到一个文件;默认输出到trimmed reads结果文件
# --max-n=COUNT：#reads中N的数量，设定整数或小数(N的占比)

# 双端测序参数
# -A ADAPTER：  #第二条reads 3'adapter
# -G ADAPTER：#第二条reads 5'adapter
# -U LENGTH： #从第二条reads上修剪的长度
# -p FILE, --paired-output=FILE： #第二条reads的输出结果
# --untrimmed-paired-output=FILE：#第一条reads没有adapter时，将第二条reads输出到文件；默认输出到trimmed reads结果文件 

# how to use this scripts
# 3.cut_adapter_PE.sh <fastq.gz file path>