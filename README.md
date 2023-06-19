# RNA-seq-2023
This is my RNA-seq  analysis update library

``` author:zhang jian ```
``` date:2023.6.19 ```
``` version:1.0v ```
``` e-mail:zhangjian199567@outlook.com ```

### this is my RNA-seq analysis protocol (SHELL + R)
### THIS LIBRARY INCLUDE 2 PART
  1.RNA-SEQ UPSTEAM ANALYSIS
    1) RAW DATA QC scripts；
    2）raw DATA QC clean scripts;
    3) clean data  QC scripts;
    4) remove MT-RNA scripts(option);
    5) mapping to genome by (HISAT2 or STAR);
    6) count by featurecount.
  
  2.RNA-SEQ DOWNSTEAM ANALYSIS 
    1) general DEG analysis R scripts；
    2）GO KEGG GSEA general R scripts;
    3) other analysis scripts
    
### need run env
#### shell env softaqre:
    1) fastqc
    2) multiqc
    3) cutadapt
    4) hisat2
    5) STAR
    6) featureCounts

#### R package
```
crayon,praise,progress,DESeq2,dbplyr,FactoMineR,plyr,stringr,Hmisc,biomaRt,ggridges,latex2exp,ggpubr,ggthemes,ggrepel,rtracklayer,ggridges,RColorBrewer,pheatmap,dplyr.plotly,htmlwidgetsshiny,reshape2,clusterProfiler,topGO,dbplyr,pathview,DOSE,org.Mm.eg.db,tidyr,dplyr,forcats,ggplot2,ggpubr,enrichplot,cowplot,stringr,RColorBrewer

```

### RNA-SEQ UPSTEAM ANALYSIS
```
# 1. downloda data form paper or seq company 
  you can use RNA-SEQ-UPSTEAM folder STEP1.2.Dowmload_data_form_EMBI-EBI.sh or STEP1.1.Dowmload_data_huawei_OBS.sh scripts 

# 2.create 1.raw_data folder & move you fq.gz file in this folder

# 3.use FASTQC & MULTIQC QC raw-data
  you can use RNA-SEQ-UPSTEAM folder STEP2.1.Data_1st_qc_by_fastqc.sh scripts
  
# 4.use cutadapter clean data
  you can use RNA-SEQ-UPSTEAM folder STEP2.2.Cutadapt_raw_data_by_cutadapt_1.0V.sh & STEP2.3.Readapter_data_qc_by_fastqc.sh remove low quality read and adapter & qc clean data result
  
# 5.use cutadapt remove 5' 15bp read (option)
 you can use RNA-SEQ-UPSTEAM folder STEP2.4.Cutadapt_option_rm_5_15bp_by_cutadapt_1.0V.sh & STEP2.5.Rm5_15bp_fastqc_report.sh
 
# 6.remove MT GENOME read if you qc result have MT-RNA-READ YOU should use down scripts remove mt-rna
  RNA-SEQ-UPSTEAM folder  STEP3.1.Build_mt-RNA-index.sh & STEP3.2.Remove_mapping_mt-RNA_hisat2.sh
  (PS:for reference genome file you can change & download for NCBI OR ENSMBL)
  
# 7.build index & mapping by HISAT2 OR STAR
  YOU can use RNA-SEQ-UPSTEAM folder include mapping or build tag scripts to mapping read to genome
  # NOTE: IF YOU BUILDING LIBRARY IS Strand-specific library preparation,you should use include dutp tag scripts
  
# 8.featurecount 
  YOU can use RNA-SEQ-UPSTEAM folder include featurecount to count featurecount 
  # NOTE: IF YOU BUILDING LIBRARY IS Strand-specific library preparation,you should use include dutp tag scripts
  # NOTE: YOU CAN USE 3 level to count RNA-seq featurecount

```

