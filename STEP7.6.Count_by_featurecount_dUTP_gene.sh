# author:zhangjian
# date:2023.1.13
# This is an count scripts for rna-seq (only suit dutp rna-seq)

# set run env
root_dir=$1
count_dir=7.count_gene_dutp
num_threads=10
annotation_file=/data/genome/new_mus/ensembl_mm39/annotation/Mus_musculus.GRCm39.109.gtf
# change dir
cd $root_dir
# create dir
mkdir -p $1/$count_dir
# use featurecount count read
for i in ./6.mapping_result/genome_mapping/*Aligned.sortedByCoord.out.bam
do
file_name=`basename $i Aligned.sortedByCoord.out.bam`
featureCounts  -F GTF \
               -t gene \
               -g gene_id \
               -T $num_threads \
               -s 2 \
               -Q 15 \
               -p \
               -a $annotation_file \
               -o ./$count_dir/${file_name}_gene.read.count.tsv \
               ./6.mapping_result/genome_mapping/${file_name}Aligned.sortedByCoord.out.bam
done
echo "count finish!!!!!!"

# use multqc to qc featurecount 
multiqc ./$count_dir -o ./$count_dir