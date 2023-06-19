# author:zhangjian
# date:2023.1.13
# This is an count scripts for rna-seq 
# set run env
root_dir=$1
count_dir=8.count_exon
num_threads=10
annotation_file=/data/genome/new_mus/ensembl_mm39/annotation/Mus_musculus.GRCm39.109.gtf

# change dir
cd $root_dir
# create dir
mkdir -p $1/$count_dir

# use featurecount count read
for i in ./7.bam_file/*.sort.bam
do
file_name=`basename $i .sort.bam`
featureCounts  -F GTF \
               -t exon \
               -g gene_id \
               -T $num_threads \
               -Q 15 \
               -p \
               -a $annotation_file \
               -o ./$count_dir/${file_name}_exon.read.count.tsv \
               ./7.bam_file/${file_name}.sort.bam
done
echo "count finish!!!!!!"

# use multqc to qc featurecount 
multiqc ./$count_dir -o ./$count_dir