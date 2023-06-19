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
crayon,praise,progress,DESeq2,dbplyr,FactoMineR,plyr,stringr,Hmisc,biomaRt,ggridges,latex2exp,ggpubr,ggthemes,ggrepel,rtracklayer,ggridges,RColorBrewer,pheatmap,dplyr.plotly,htmlwidgetsshiny,reshape2,clusterProfiler,topGO,dbplyr,pathview,DOSE,org.Mm.eg.db,tidyr,dplyr,forcats,ggplot2,ggpubr,enrichplot,cowplot,stringr,RColorBrewer
