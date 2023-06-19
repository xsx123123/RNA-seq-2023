#!/usr/bin/env R
# author:zhang jian
# date:2023.6.16
# version:5.2.0v
# e-mail:zhangjian199567@outlook.com
# This is an rna-seq downsteam analysis scripts for rat
##########################################################################
rm(list=ls())
################################# PATH-1 #################################
# PATH-1 pre define function & pre loading package
suppressMessages(library("crayon"))
suppressMessages(library("praise"))
suppressMessages(library("progress"))
# Function1:print run condition style 1
print_color_note_type1 <- function(sologo){
  cat(underline(bold(bgCyan("#-----------------------------------------------#####-----------------------------------------------#"))))
  cat("\n")
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                                           "),sologo,"\n")
  cat(rep(" ",17),factor)
  cat("\n")
  cat("\n")
}
# Function2:print run condition style 2
print_color_note_type2 <- function(sologo){
  cat("\n")
  cat("\n")
  cat(underline(bold(bgCyan("#-----------------------------------------------#####-----------------------------------------------#"))))
  cat("\n")
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                         "),sologo,"\n")
  cat(rep(" ",18),factor)
  cat("\n")
  cat(underline(bold(bgCyan("#-----------------------------------------------#####-----------------------------------------------#"))))
  cat("\n")
  cat("\n")
}
# Function2:print run condition style 3
print_color_note_type2_warring <- function(sologo){
  cat("\n")
  cat("\n")
  cat(underline(bold(bgRed("#-----------------------------------------------#####-----------------------------------------------#"))))
  cat("\n")
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                "),sologo,"\n")
  cat(rep(" ",8),factor)
  cat("\n")
  cat(underline(bold(bgRed("#-----------------------------------------------#####-----------------------------------------------#"))))
  cat("\n")
  cat("\n")
}
# Function3:print run condition style 4
print_color_note_type3 <- function(sologo){
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                                           "),sologo,"\n")
  cat(rep(" ",18),factor)
  cat("\n")
  cat(underline(bold(bgCyan("#-----------------------------------------------#####-----------------------------------------------#"))))
  cat("\n")
  cat("\n")
}
# Function4:print ProgressBar
ProgressBar <- function(){
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total = 6, clear = FALSE, width = 120
  )
  
  for (i in 1:6) {
    pb$tick()
    Sys.sleep(0.5)
  }
}
# Function5:print run parameter
get_run_parameter <- function(){
  run_parameter <- c("DEG analysis root_dir is :","featurecount count dir is :","sample infor file of DEG is :",
                     "special file dir is :","DEG result save dir is :","ref (backgroud) sample group name is:",
                     "featurecount count level is :","special gene id is :","paper Volcano plot group infor:","paper Volcano plot save name&title:")
  run_parameter <- as.data.frame(run_parameter)
  run_parameter$parameter <- c(root_dir,count_dir,
                               sample_file_name_dir,
                               special,output_name,
                               ref_name,count_level,
                               gene_1,EXP_NAEE,figure_name)
  print_color_note_type1("place check run parameter!!!")
  print(run_parameter)
  print_color_note_type3("place check run parameter!!!")
  cat("\n")
  ProgressBar()
  cat("\n")
}
# Function6:desktop ensmbl path
desktop_get_ensmbl_path <- function(){
  if (special=="rat"){
    ensmbl_path <- "/Users/zhangjian/Desktop/Pre-doctoral/local_database/ensmbl/rat_ensembl_gene_id_mgi_symbol_2023-05-08_11:05:42.csv" 
  }else{
    if(special=="mus"){
      ensmbl_path <- "/Users/zhangjian/Desktop/Pre-doctoral/local_database/ensmbl/mus_ensembl_gene_id_mgi_symbol_2023-05-07_16:49:09.csv" 
    }else{
      if(special=="human"){
        ensmbl_path <- "/Users/zhangjian/Desktop/Pre-doctoral/local_database/ensmbl/human_ensembl_gene_id_mgi_symbol_2023-05-07_16:49:09.csv"
      }
    }
  }
  return(ensmbl_path)
}
# Function7:severs ensmbl path
severs_get_ensmbl_path <- function(){
  if (special=="rat"){
    ensmbl_path <- "/data/genome/ensmbl_genome_annotation_online/rat_ensembl_gene_id_mgi_symbol_2023-05-08_11:05:42.csv" 
  }else{
    if(special=="mus"){
      ensmbl_path <- "/data/genome/ensmbl_genome_annotation_online/mus_ensembl_gene_id_mgi_symbol_2023-04-21_20:47:17.csv" 
    }else{
      if(special=="human"){
        ensmbl_path <- "/data/genome/ensmbl_genome_annotation_online/human_ensembl_gene_id_mgi_symbol_2023-05-07_16:49:09.csv"
      }
    }
  }
  return(ensmbl_path)
}
# Function8:auto get ensmbl path
desktop_severs_ensmbl <- function(){
  path <- getwd()
  tag <- strsplit(path,"/")[[1]][2]
  if (tag == "Users"){
    ensmbl_path <- desktop_get_ensmbl_path()
  }else{
    if (tag=="data"){
      ensmbl_path <-severs_get_ensmbl_path()
    }
  }
  return(ensmbl_path)
}
# Function9:check package install
install_package_check <- function(){
  list <- as.data.frame(installed.packages())
  # set need install package
  need_package_list <- c("crayon","praise","progress",
                         "DESeq2","dbplyr","FactoMineR",
                         "plyr","stringr","Hmisc",
                         "biomaRt","ggridges","latex2exp",
                         "ggpubr","ggthemes","ggrepel",
                         "rtracklayer","ggridges","RColorBrewer",
                         "pheatmap","dplyr","plotly",
                         "htmlwidgets","shiny","reshape2")
  for (i in need_package_list){
    if(i %in% list$Package){
      cat("Analysis package",i,"install !!!!   ")
      cat(bold(green("PASS!!!")))
      cat("\n")
    }else{
      cat("Analysis package",i,"not install !!!!")
      cat(bold(red("FAIL!!!")))
      stop("you should install depend package!!! !!!!")
    }
    Sys.sleep(0.5)
  }
}
################################# PATH-2 #################################
# PATH-2 loading shell command
# print run condition
print_color_note_type1("loading shell command DO!!!")
args<-commandArgs(TRUE)
# deg aanalysis root dir
root_dir <- args[1]
# count file dir by featurecount
count_dir <- args[2]
# group infor for DEG
sample_file_name_dir <- args[3]
# set special name:rat mus human
special <- "mus"
# get ensmbl path
ensmbl_path <- desktop_severs_ensmbl()
# DEG  result folder
output_name <- "8.deg_exon"
# DEG reference group NAME
ref_name <- "WT"
# count level by featurecount
count_level <- "exon"
# need add special tag
gene_1 <- "Mtco3"
# set group infor for  Volcano plot  of paper
EXP_NAEE <- "Mut MT-CO3"
# Volcano plot  of paper save name and title
figure_name <- "C6 cell Mut MT-CO3 vs WT "
# print run parameter
get_run_parameter()
# print run condition
print_color_note_type3("loading shell command done DONE!!!")
################################# PATH-3 #################################
# PATH-3 loading package
# print run condition
print_color_note_type1("loading package DO!!!")
install_package_check()
# loading package
suppressMessages(library("DESeq2"))
suppressMessages(library("dbplyr"))
suppressMessages(library("FactoMineR"))
suppressMessages(library("plyr"))
suppressMessages(library("stringr"))
suppressMessages(library("Hmisc"))
suppressMessages(library("biomaRt"))
suppressMessages(library("ggridges"))
suppressMessages(library("latex2exp"))
suppressMessages(library("ggpubr"))
suppressMessages(library("ggthemes"))
suppressMessages(library("ggrepel"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("ggridges"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("pheatmap"))
suppressMessages(library("dplyr"))
suppressMessages(library("plotly"))
suppressMessages(library("htmlwidgets"))
suppressMessages(library("shiny"))
suppressMessages(library("reshape2"))
# print run condition
print_color_note_type3("loading package done DONE!!!")
################################# PATH-4 #################################
# PATH-4 define function
# print run condition
print_color_note_type1("loading function DO!!!")
# Function1:define read featurecount file function
read_feature <- function(featurecount_path,name_file){
  name_ID <- strsplit(name_file,split='_exon.read.count.tsv')[[1]][1]
  count_tmp <- read.csv(file.path(count_dir,name_file),header = F,sep = "\t")
  count_tmp <- count_tmp[-c(1:2),]
  count_tmp <- count_tmp[,-c(2:6)]
  colnames(count_tmp) <- c("gene_id",name_ID)
  #head(count_tmp)
  return(count_tmp)
}
# Function2:define merge function for count file for featurecount
merge_count <- function(count_dir,sample_name){
  df <- read_feature(count_dir,sample_name[1])
  for (i in sample_name){
    count_tmp <- read_feature(count_dir,i)
    df <- left_join(df,count_tmp,by="gene_id")
  }
  df <- df[,-2]
  return(df)
}
# Function3:define DEG main function model
DEG <- function(count,infor,ref,dds_file_name) {
  # input data to deseq2
  infor$condition <- factor(infor$condition)
  infor$patch <- factor(infor$patch)
  infor$lane <- factor(infor$lane)
  # count <- count+1
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = infor,
                                design = ~ condition)
  # set ref group
  dds$condition <- relevel(dds$condition, ref = ref)
  # DEG
  dds <- DESeq(dds)
  # result
  res <- results(dds)
  # summary resulit
  summary(res)
  # extert DEG result
  result_data <- merge(as.data.frame(res),
                       as.data.frame(counts(dds,normalized=TRUE)),
                       by='row.names',
                       sort=FALSE)
  # GENE ID CONVERT
  # get pre convert gene id
  gene_list <- result_data$Row.names
  # gene id convert
  convert_id_list <- convert_name(gene_list)
  # rename colnames
  colnames(convert_id_list) <- c("Row.names","SYMBOL")
  # remove dup name
  convert_id_list_re <- convert_id_list[!duplicated(convert_id_list$Row.names),]
  # add symbol to DEG result
  result_data_add_name <- left_join(result_data,convert_id_list_re,by="Row.names")
  # NA data deal
  result_data_add_name <- merge_gene_id(result_data_add_name)
  # save dds result
  write.csv(result_data_add_name,file.path(output_dir,dds_file_name))
  # print dds file
  cat("DEG result file path:",file.path(output_dir,dds_file_name))
  return(result_data_add_name) # return 为将结果从函数中导出
}
# Function4:get parent folder name
get_parent_name <- function(){
  now_dir <- getwd()
  temp_list <- strsplit(now_dir,split='/')[[1]][-1]
  name <- temp_list[length(temp_list)]
  return(name)
}
# Function5:create dir
create_dir <- function(list_dir){
  for (i in list_dir) {
    if(! dir.exists(i)){
      dir.create(i)
    }
  }
  cat("create dir finish!!!!")
}
# Function6:gene id convert
convert_name <- function(gene_list){
  gene_list <- as.data.frame(gene_list);colnames(gene_list) <- "ensembl_gene_id"
  head(gene_list)
  dim(gene_list)
  rat_ensmbl <- read.csv(ensmbl_path,row.names=1)
  head(rat_ensmbl)
  merge_data <- left_join(gene_list,rat_ensmbl)
  dim(merge_data)
  # remove dup
  dup <- duplicated(merge_data$ensembl_gene_id)
  merge_data <- merge_data[!dup, ]
  dim(merge_data)
  # table NA data
  # library(org.Mm.eg.db)
  # columns(org.Mm.eg.db)
  # head(keys(org.Mm.eg.db,"ENSEMBL"))
  # head(keys(org.Mm.eg.db,"SYMBOL"))
  # id_converted_list = bitr(gene_list, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Mm.eg.db)
  return(merge_data)
}
# Function7:merge gene id
merge_gene_id <- function(result_data){
  list <- which(is.na(result_data$SYMBOL))
  for (i in list){
    result_data[as.numeric(i),]$SYMBOL <- result_data[as.numeric(i),]$Row.names
  }
  result_data_merge_id <- result_data
  return(result_data_merge_id)
}
# Function8:get up gene list
get_up_gene_padj <- function(deg_result){
  deg_result_up <- subset(deg_result,deg_result$log2FoldChange > 1 &deg_result$padj < 0.05)
  write.csv(deg_result_up,file.path(output_dir,paste0(get_parent_name(),"_DEG_padj_up_output_result.csv")))
  return(deg_result_up)
}
get_up_gene_pvalueue <- function(deg_result){
  deg_result_up <- subset(deg_result,deg_result$log2FoldChange > 1 & deg_result$pvalue < 0.05)
  write.csv(deg_result_up,file.path(output_dir,paste0(get_parent_name(),"_DEG_pvalueue_up_output_result.csv")))
  return(deg_result_up)
}
# Function9:get down gene list
get_down_gene_padj <- function(deg_result){
  deg_result_down <- subset(deg_result,deg_result$log2FoldChange < -1 & deg_result$padj < 0.05)
  write.csv(deg_result_down,file.path(output_dir,paste0(get_parent_name(),"_DEG_padj_down_output_result.csv")))
  return(deg_result_down)
}
get_down_gene_pvalueue <- function(deg_result){
  deg_result_up <- subset(deg_result,deg_result$log2FoldChange < -1 &deg_result$pvalue < 0.05)
  write.csv(deg_result_up,file.path(output_dir,paste0(get_parent_name(),"_DEG_pvalueue_down_output_result.csv")))
  return(deg_result_up)
}
# Function10:cut sample name
cut_name <- function(pcaData){
  name_list <- pcaData$name
  cut_name_list <- c()
  for (i in name_list){
    name <- strsplit(i,split='_L[1-9]_')[[1]][1]
    cut_name_list <- append(cut_name_list,name)
  }
  return(cut_name_list )
}
# Function11:draw PCA  plot
plot_PCA <- function(count,infor,dds_file_name) {
  # dds_file_name <- paste0(get_parent_name(),"_PCA_plot")
  infor <- sample_infor
  infor$condition <- factor(infor$condition)
  # count <- count+1
  dds <- DESeqDataSetFromMatrix(countData = df_1,colData = sample_infor,
                                design = ~ condition)
  vsd <- vst(dds, blind=FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
  cut_name_list <- cut_name(pcaData)
  pcaData$name_2 <- cut_name_list
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  axis_range <- max(abs(pcaData$PC1),abs(pcaData$PC2))
  p <- ggplot(data = pcaData, aes(x = PC1,y = PC2)) +  
    geom_point(aes(fill = condition),size=3,stroke = 0.6,alpha=0.6,shape=21) +      
    scale_fill_manual(values = c('#ffb6b9','#8ac6d1',"#9dd3a8")) +  
    geom_vline(xintercept=0,lty=4,
               col="black",lwd=0.15) +          
    geom_hline(yintercept = 0,lty=4,
               col="black",lwd=0.15) +          
    geom_text_repel(aes(label = name_2), size = 3, show.legend = FALSE, box.padding = unit(0.8,'lines'),
                    size = 1.5,segment.alpha = 0.4,segment.size = 0.3,segment.color = "black") +    
    scale_x_continuous(limits=c(-(axis_range),axis_range),n.breaks = 8) +           
    scale_y_continuous(limits=c(-(axis_range),axis_range),n.breaks = 8) +          
    theme(panel.grid = element_blank(),legend.position = 'bottom',,panel.background = element_rect(color = 'black',fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent')) +  
    labs(x =  paste('PCA1:',percentVar[1], '%'),
         y = paste('PCA2:',percentVar[2], '%'),color = '')  + guides(fill=guide_legend(nrow=1,title=NULL,override.aes = list(size=3,alpha=0.5),direction="horizontal"))
  
  ggsave(paste0(dds_file_name,".pdf"),plot = p,device = "pdf",path = figure_dir,width = 10,height = 11,units = "cm")
  ggsave(paste0(dds_file_name,".png"),plot = p,device = "png",path = figure_dir,width = 10,height = 11,units = "cm",dpi=1000)
  cat("PCA plot save path:",file.path(figure_dir,dds_file_name))
}
# Function12:cut sample name-2
cut_name_2 <- function(deg_result){
  list <- colnames(deg_result)
  cut_name_list <- c()
  for (i in list){
    name <- strsplit(i,split='_L[1-9]_')[[1]][1]
    cut_name_list <- append(cut_name_list,name)
  }
  return(cut_name_list )
}
# function13:cut sampleinfor file
cut_name_3 <- function(sample_infor){
  list <- sample_infor$group_name
  cut_name_list <- c()
  for (i in list){
    name <- strsplit(i,split='_L[1-9]_')[[1]][1]
    cut_name_list <- append(cut_name_list,name)
  }
  
  sample_infor$condition <- cut_name_list
  return(sample_infor)
}
# Function14:remove sample point
remove_sample_point <- function(deal_name){
  for (i in c(1:length(deal_name))){
    if (grepl(".", deal_name[i], fixed=TRUE)){
      deal_name[i] <- gsub("\\.", "-", deal_name[i])
    }else{
      deal_name[i] <- deal_name[i]
    }
  }
  return(deal_name)
}
# Function15:plot heatmap cor
plot_heatmap_corr <- function(deg_result){
  sample_infor <- cut_name_3(sample_infor)
  data <- deg_result
  rownames(data) <- data$Row.names
  data <- data[,-c(1:7)]
  data <- data[,-(nrow(sample_infor)+1)]
  name <- cut_name_2(data)
  colnames(data) <- name
  res2_spearman <- rcorr(as.matrix(data),type = "spearman")
  cor_maritx <- res2_spearman$r
  cor_maritx <- data.frame(cor_maritx)
  colnames(cor_maritx) <- rownames(cor_maritx)
  cor_maritx$id <- rownames(cor_maritx)
  #convert long data frame
  temp <- data.frame()
  for (i in c(1:dim(cor_maritx)[1])){
    # get row
    temp_row <- cor_maritx[,c(i,dim(cor_maritx)[2])]
    temp_row$name <- colnames(temp_row)[1]
    colnames(temp_row) <- c("cor","group_1","group_2")
    temp <- rbind(temp,temp_row)
  }
  # remove saame group name
  list <- !duplicated(sample_infor$condition)
  remove_dup_list <- sample_infor$condition[list]
  # sample name deal
  temp$group_1 <- remove_sample_point(temp$group_1)
  temp$group_2 <- remove_sample_point(temp$group_2)
  # set factor for everyone sample name
  temp$group_2 <- factor(temp$group_2,levels=remove_dup_list)
  temp$group_1 <- factor(temp$group_1,levels=remove_dup_list)
  p <- ggplot(temp,aes(x=group_1,y=group_2))+ 
    geom_tile(aes(fill=cor))+
    scale_fill_gradient(low="#B1E3FA",high = "#2980B9",name=NULL)+
    labs(x="",y="")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "right",
          legend.title=element_text(size = 10), 
          text = element_text(size = 9,family="sans"),title = element_text(size = 7),           
          axis.title.x = element_text(size = 7),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 8,angle=45,hjust=0.9,vjust=1),axis.text.y = element_text(size = 8),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1),
          panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) 
  ggsave(paste0(get_parent_name(),"_cor_value.pdf"),plot=p,width = 12,height = 10,units = "cm",path = figure_dir)
  ggsave(paste0(get_parent_name(),"_cor_value.png"),plot=p,device = "png",width = 12,height = 10,units = "cm",path = figure_dir,dpi=1000)
  relative_value <- res2_spearman$r
  write.csv(relative_value,file.path(output_dir,paste0(get_parent_name(),"_cor_heatmapvalueue.csv")))
  relative_pvalueue <- res2_spearman$P
  write.csv(relative_pvalueue,file.path(output_dir,paste0(get_parent_name(),"_cor_heatmap_pvalueue.csv")))
  cat("\n")
  cat("heatmap cor plot save path:",file.path(output_dir,paste0(get_parent_name(),"_cor_pvalueue.csv")))
  cat("\n")
}
# Function16:plot Volcano Plot
draw_volcano <- function(deg_result){
  deg_result <- deg_result
  deg_result$log10 <- -log10(deg_result$pvalue)
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$pvalue < 0.05) & (deg_result$log2FoldChange > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$pvalue  < 0.05) & (deg_result$log2FoldChange < -1))] = "Down-regulated"
  deg_result <- deg_result[order(deg_result$pvalue),]
  up.genes <- head(deg_result$gene_name[which(deg_result$Group == "Up-regulated")],10)
  down.genes <- head(deg_result$gene_name[which(deg_result$Group == "Down-regulated")],10)
  deg_result.top.genes <- c(as.character(up.genes), as.character(down.genes))
  deg_result$label[match(deg_result.top.genes,deg_result$gene_name)] <- deg_result.top.genes
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),10)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),10)
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  p <- ggplot(deg_result, aes(x = log2FoldChange, y = log10, colour=Group)) +  
    geom_point(alpha=1, size=1.3,alpha=0.7) +
    geom_point(data = deg_result_up,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    geom_point(data = deg_result_down,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    scale_color_manual(values=c("#56959A","#C7C7C7","#9C56CD")) + 
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.15) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.15) +  
    labs(x="Fold Change (log2)",y= TeX(r"(\textit{P} value (-log10) )"),title =paste0(get_parent_name()," Volcano Plot")) + 
    geom_text_repel(data = deg_result_up,aes(log2FoldChange, log10, label= SYMBOL),size=2,colour="black",
                    segment.alpha = 0.5,segment.size = 0.3,segment.color = "black",
                    box.padding=unit(0.5, "lines"),point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(0.5, 6),max.overlaps = 25) + 
    geom_text_repel(data = deg_result_down,aes(log2FoldChange, log10, label= SYMBOL),size=2,colour="black",
                    segment.alpha =0.5,segment.size = 0.3,segment.color = "black",
                    box.padding=unit(0.5, "lines"),point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(-6, -0.5),max.overlaps = 25) + 
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),
                       n.breaks = 10) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 8) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),legend.position="none",
          legend.title=element_text(size = 10,face = "normal",family="Times"), 
          text = element_text(size = 9,family="sans"),title = element_text(size = 7),           
          axis.title.x = element_text(size = 7),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))        
  # 1.5 save Volcano plot
  ggsave(file.path(figure_dir,paste0(get_parent_name()," Volcano Plot.pdf")),plot = p,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," Volcano Plot.png")),device = "png",plot = p,width = 9,height = 10,units = "cm",dpi = 1000)
  cat("volcan figure save path:",file.path(figure_dir,paste0(get_parent_name(),"Volcano Plot.pdf/.png")))
  return(deg_result)
}
# Function17:plot Volcano Plot
draw_volcano_padj <- function(deg_result){
  deg_result <- deg_result
  deg_result$log10 <- -log10(deg_result$padj)
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$padj < 0.05) & (deg_result$log2FoldChange > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$padj  < 0.05) & (deg_result$log2FoldChange < -1))] = "Down-regulated"
  deg_result <- deg_result[order(deg_result$padj),]
  up.genes <- head(deg_result$gene_name[which(deg_result$Group == "Up-regulated")],10)
  down.genes <- head(deg_result$gene_name[which(deg_result$Group == "Down-regulated")],10)
  deg_result.top.genes <- c(as.character(up.genes), as.character(down.genes))
  deg_result$label[match(deg_result.top.genes,deg_result$gene_name)] <- deg_result.top.genes
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),10)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),10)
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  p <- ggplot(deg_result, aes(x = log2FoldChange, y = log10, colour=Group)) +  
    geom_point(alpha=1, size=1.3,alpha=0.7) +
    geom_point(data = deg_result_up,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    geom_point(data = deg_result_down,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    scale_color_manual(values=c("#56959A","#C7C7C7","#9C56CD")) + 
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.15) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.15) +  
    labs(x="Fold Change (log2)",y= TeX(r"(\textit{P} adj (-log10) )"),title =paste0(get_parent_name()," Volcano Plot")) + 
    geom_text_repel(data = deg_result_up,aes(log2FoldChange, log10, label= SYMBOL),size=2,colour="black",
                    segment.alpha = 0.5,segment.size = 0.3,segment.color = "black",
                    box.padding=unit(0.5, "lines"),point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(0.5, 6),max.overlaps = 25) + 
    geom_text_repel(data = deg_result_down,aes(log2FoldChange, log10, label= SYMBOL),size=2,colour="black",
                    segment.alpha =0.5,segment.size = 0.3,segment.color = "black",
                    box.padding=unit(0.5, "lines"),point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(-6, -0.5),max.overlaps = 25) + 
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),
                       n.breaks = 10) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 8) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),legend.position="none",
          legend.title=element_text(size = 10,face = "normal",family="Times"), 
          text = element_text(size = 9,family="sans"),title = element_text(size = 7),           
          axis.title.x = element_text(size = 7),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))        
  # 1.5 save Volcano plot
  ggsave(file.path(figure_dir,paste0(get_parent_name()," Volcano Plot.pdf")),plot = p,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," Volcano Plot.png")),device = "png",plot = p,width = 9,height = 10,units = "cm",dpi = 1000)
  cat("volcan figure save path:",file.path(figure_dir,paste0(get_parent_name(),"Volcano Plot.pdf/.png")))
  return(deg_result)
}
# Function17: draw volcano style 2
draw_volcano_specail <- function(deg_result,gene_1,plot_name,dir){
  data <- deg_result
  data$FC <- 2^(data$log2FoldChange)
  data$log10 <- -log10(data$pvalue)
  pvalue = 0.05
  log2FC = 1
  f <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_1 <- data.frame(x = inputx + log2FC, y = y)
    return(dff_1)
  }
  r <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_2 <- data.frame(x = -(inputx + log2FC), y = y)
    return(dff_2)
  }
  line_data_1 <- f(5)
  line_data_2 <- r(5)
  data$curve_y <- case_when(
    data$log2FoldChange > 0 ~ 1/(data$log2FoldChange-log2FC) + (-log10(pvalue)),
    data$log2FoldChange <= 0 ~ 1/(-data$log2FoldChange-log2FC) + (-log10(pvalue)))
  data$group2 <- case_when(
    data$log10 > data$curve_y & data$log2FoldChange >= log2FC ~ 'Up-regulated',
    data$log10 > data$curve_y & data$log2FoldChange <= -log2FC ~ 'Down-regulated',
    TRUE ~ 'Non-regulated')
  # head(data)
  GENE <-subset(data,data$SYMBOL==gene_1)
  data_up_gene <- subset(data,data$group2=="Up-regulated")
  data_up_gene_top10 <- head(data_up_gene[order(data_up_gene$pvalue),],10)
  data_down_gene <- subset(data,data$group2=="Down-regulated")
  data_down_gene_top10 <- head(data_down_gene[order(data_down_gene$pvalue),],10)
  # get y_aes
  y_aes <- data$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- data$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw plot 
  p <- ggplot(data,aes(x=log2FoldChange,y=log10,color=group2))+
    geom_point(size=1.5,stroke = 0.5,,alpha=0.2)+ 
    geom_point(data = GENE,aes(log2FoldChange, log10,),alpha=0.8, size=2.4,shape = 1,color="black")+
    theme_linedraw()+
    geom_line(data = line_data_1,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    geom_line(data = line_data_2,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    scale_color_manual(values=c("#56959A","#C7C7C7","#9C56CD")) + 
    labs(x="Fold Change (log2)",y= TeX(r"(\textit{P} value (-log10) )"),title =paste0(get_parent_name()," Volcano Plot")) + 
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    geom_text_repel(data = GENE,aes(log2FoldChange, log10, label= SYMBOL),size=1.5,colour="black",
                    segment.alpha = 0.5,segment.size = 0.3,segment.color = "black",min.segment.length = 0,
                    box.padding = 0.5,point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,,max.overlaps = 25) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),
                       n.breaks = 10) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 8) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),legend.position="none",
          legend.title=element_text(size = 10,face = "normal",family="Times"), 
          text = element_text(size = 9,family="sans"),title = element_text(size = 7),           
          axis.title.x = element_text(size = 7),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))
  # 1.5 save Volcano plot
  ggsave(file.path(dir,paste0(get_parent_name()," ",gene_1," ",plot_name," Volcano Plot special.pdf")),plot = p,width = 9,height = 10,units = "cm")
  ggsave(file.path(dir,paste0(get_parent_name()," ",gene_1," ",plot_name," Volcano Plot special.png")),device = "png",plot = p,width = 9,height = 10,units = "cm",dpi = 1000)
  # cat("\n")
  # cat("    ","\n","🍻🍻 (●´∀｀●)ﾉ draw ",gene_1,"FINISH (●´∀｀●)ﾉ 🍻🍻!!!!")
  # cat("\n")
  # cat("volcan figure save path:",file.path(dir,paste0(get_parent_name()," ",gene_1," ",plot_name," Volcano Plot special.pdf/png")))
  # cat("\n")
}
# Function18: draw volcano style 2
draw_volcano_specail_padj <- function(deg_result,gene_1,plot_name,dir){
  data <- deg_result
  data$FC <- 2^(data$log2FoldChange)
  data$log10 <- -log10(data$padj)
  pvalue = 0.05
  log2FC = 1
  f <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_1 <- data.frame(x = inputx + log2FC, y = y)
    return(dff_1)
  }
  r <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_2 <- data.frame(x = -(inputx + log2FC), y = y)
    return(dff_2)
  }
  line_data_1 <- f(5)
  line_data_2 <- r(5)
  data$curve_y <- case_when(
    data$log2FoldChange > 0 ~ 1/(data$log2FoldChange-log2FC) + (-log10(pvalue)),
    data$log2FoldChange <= 0 ~ 1/(-data$log2FoldChange-log2FC) + (-log10(pvalue)))
  data$group2 <- case_when(
    data$log10 > data$curve_y & data$log2FoldChange >= log2FC ~ 'Up-regulated',
    data$log10 > data$curve_y & data$log2FoldChange <= -log2FC ~ 'Down-regulated',
    TRUE ~ 'Non-regulated')
  # head(data)
  GENE <-subset(data,data$SYMBOL==gene_1)
  data_up_gene <- subset(data,data$group2=="Up-regulated")
  data_up_gene_top10 <- head(data_up_gene[order(data_up_gene$padj),],10)
  data_down_gene <- subset(data,data$group2=="Down-regulated")
  data_down_gene_top10 <- head(data_down_gene[order(data_down_gene$padj),],10)
  # get y_aes
  y_aes <- data$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- data$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw plot 
  p<-ggplot(data,aes(x=log2FoldChange,y=log10,color=group2))+
    geom_point(size=1.5,stroke = 0.5,,alpha=0.2)+ 
    geom_point(data = GENE,aes(log2FoldChange, log10,),alpha=0.8, size=2.4,shape = 1,color="black")+
    theme_linedraw()+
    geom_line(data = line_data_1,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    geom_line(data = line_data_2,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    scale_color_manual(values=c("#56959A","#C7C7C7","#9C56CD")) + 
    labs(x="Fold Change (log2)",y= TeX(r"(\textit{P} adj (-log10) )"),title =paste0(get_parent_name()," Volcano Plot")) + 
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    geom_text_repel(data = GENE,aes(log2FoldChange, log10, label= SYMBOL),size=1.5,colour="black",
                    segment.alpha = 0.5,segment.size = 0.3,segment.color = "black",min.segment.length = 0,
                    box.padding = 0.5,point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,,max.overlaps = 25) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),
                       n.breaks = 10) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 8) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),legend.position="none",
          legend.title=element_text(size = 10,face = "normal",family="Times"), 
          text = element_text(size = 9,family="sans"),title = element_text(size = 7),           
          axis.title.x = element_text(size = 7),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))
  # 1.5 save Volcano plot
  ggsave(file.path(dir,paste0(get_parent_name()," ",gene_1," ",plot_name," Volcano Plot special.pdf")),plot = p,width = 9,height = 10,units = "cm")
  ggsave(file.path(dir,paste0(get_parent_name()," ",gene_1," ",plot_name," Volcano Plot special.png")),device = "png",plot = p,width = 9,height = 10,units = "cm",dpi = 1000)
  # cat("\n")
  # cat("    ","\n","🍻🍻 (●´∀｀●)ﾉ draw ",gene_1,"FINISH (●´∀｀●)ﾉ 🍻🍻!!!!")
  # cat("\n")
  # cat("volcan figure save path:",file.path(dir,paste0(get_parent_name()," ",gene_1," ",plot_name," Volcano Plot special.pdf/png")))
  # cat("\n")
}
# Function19: draw function
draw_volcano_up_down <- function(data,data_up_gene,data_down_gene,line_data_1,line_data_2){
  # x_aes deal
  x_aes <- data$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # get y_aes
  y_aes <- data$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # draw plot
  p <-ggplot(data,aes(x=log2FoldChange,y=log10,color=group2))+
    geom_point(size=1.5,stroke = 0.5,alpha=0.8)+ 
    geom_point(data = data_up_gene,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    geom_point(data = data_down_gene,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    theme_linedraw()+
    geom_line(data = line_data_1,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    geom_line(data = line_data_2,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    scale_color_manual(values=c("#56959A","#C7C7C7","#9C56CD")) + 
    labs(x="Fold Change (log2)",y= TeX(r"(\textit{P} value (-log10) )"),title =paste0(get_parent_name()," Volcano Plot")) + 
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    geom_text_repel(data = data_up_gene,aes(log2FoldChange, log10, label= SYMBOL),size=1.5,colour="black",
                    segment.alpha = 0.5,segment.size = 0.3,segment.color = "black",min.segment.length = 0,
                    box.padding = 0.5,point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(0.5, max(abs(data$log2FoldChange))*1.4),max.overlaps = 30) + 
    geom_text_repel(data = data_down_gene,aes(log2FoldChange, log10, label= SYMBOL),size=1.5,colour="black",
                    segment.alpha =0.5,segment.size = 0.3,segment.color = "black",min.segment.length = 0,
                    box.padding = 0.5,point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(-(max(abs(data$log2FoldChange))*1.4), -0.5),max.overlaps = 30) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),
                       n.breaks = 10) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 8) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),legend.position="none",
          legend.title=element_text(size = 10,face = "normal",family="Times"), 
          text = element_text(size = 9,family="sans"),title = element_text(size = 7),           
          axis.title.x = element_text(size = 7),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))
  return(p)  
}
# Function20: draw function
draw_volcano_up_down_padj <- function(data,data_up_gene,data_down_gene,line_data_1,line_data_2){
  # x_aes deal
  x_aes <- data$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # get y_aes
  y_aes <- data$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # draw plot
  p <-ggplot(data,aes(x=log2FoldChange,y=log10,color=group2))+
    geom_point(size=1.5,stroke = 0.5,alpha=0.8)+ 
    geom_point(data = data_up_gene,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    geom_point(data = data_down_gene,aes(log2FoldChange, log10,),alpha=0.3, size=2.4,shape = 1,color="black") +
    theme_linedraw()+
    geom_line(data = line_data_1,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    geom_line(data = line_data_2,aes(x = x, y = y),color = "black",lty = "dashed",size = 0.3,lineend = "round",alpha=0.4) +
    scale_color_manual(values=c("#56959A","#C7C7C7","#9C56CD")) + 
    labs(x="Fold Change (log2)",y= TeX(r"(\textit{P} adj (-log10) )"),title =paste0(get_parent_name()," Volcano Plot")) + 
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    geom_text_repel(data = data_up_gene,aes(log2FoldChange, log10, label= SYMBOL),size=1.5,colour="black",
                    segment.alpha = 0.5,segment.size = 0.3,segment.color = "black",min.segment.length = 0,
                    box.padding = 0.5,point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(0.5, max(abs(data$log2FoldChange))*1.4),max.overlaps = 30) + 
    geom_text_repel(data = data_down_gene,aes(log2FoldChange, log10, label= SYMBOL),size=1.5,colour="black",
                    segment.alpha =0.5,segment.size = 0.3,segment.color = "black",min.segment.length = 0,
                    box.padding = 0.5,point.padding=unit(0.3, "lines"),force = 2,max.iter = 3e3,
                    xlim=c(-(max(abs(data$log2FoldChange))*1.4), -0.5),max.overlaps = 30) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),
                       n.breaks = 10) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 8) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),legend.position="none",
          legend.title=element_text(size = 10,face = "normal",family="Times"), 
          text = element_text(size = 9,family="sans"),title = element_text(size = 7),           
          axis.title.x = element_text(size = 7),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))
  return(p)  
}
# Function21: draw volcano style 2
draw_volcano_type_2 <- function(deg_result_add_lable){
  data <- deg_result
  data$FC <- 2^(data$log2FoldChange)
  data$log10 <- -log10(data$pvalue)
  pvalue = 0.05
  log2FC = 1
  f <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_1 <- data.frame(x = inputx + log2FC, y = y)
    return(dff_1)
  }
  r <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_2 <- data.frame(x = -(inputx + log2FC), y = y)
    return(dff_2)
  }
  line_data_1 <- f(5)
  line_data_2 <- r(5)
  #新增曲线数值列：
  ##每列log2FoldChange值在曲线上对应的y轴坐标；
  data$curve_y <- case_when(
    data$log2FoldChange > 0 ~ 1/(data$log2FoldChange-log2FC) + (-log10(pvalue)),
    data$log2FoldChange <= 0 ~ 1/(-data$log2FoldChange-log2FC) + (-log10(pvalue)))
  #根据曲线新增上下调分组标签：
  data$group2 <- case_when(
    data$log10 > data$curve_y & data$log2FoldChange >= log2FC ~ 'Up-regulated',
    data$log10 > data$curve_y & data$log2FoldChange <= -log2FC ~ 'Down-regulated',
    TRUE ~ 'Non-regulated')
  head(data)
  data_up_gene <- subset(data,data$group2=="Up-regulated")
  data_up_gene_top10 <- head(data_up_gene[order(data_up_gene$pvalue),],10)
  data_up_gene_top30 <- head(data_up_gene[order(data_up_gene$pvalue),],30)
  data_down_gene <- subset(data,data$group2=="Down-regulated")
  data_down_gene_top10 <- head(data_down_gene[order(data_down_gene$pvalue),],10)
  data_down_gene_top30 <- head(data_down_gene[order(data_down_gene$pvalue),],30)
  # 1.5 save Volcano plot
  p1 <- draw_volcano_up_down(data,data_up_gene_top10,data_down_gene_top10,line_data_1,line_data_2)
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top10 Volcano Plot type_2.pdf")),plot = p1,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top10 Volcano Plot type_2.png")),device = "png",plot = p1,width = 9,height = 10,units = "cm",dpi = 1000)
  #
  p2 <- draw_volcano_up_down(data,data_up_gene_top30,data_down_gene_top30,line_data_1,line_data_2)
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(p-value).pdf")),plot = p2,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(p-value).png")),device = "png",plot = p2,width = 9,height = 10,units = "cm",dpi = 1000)
  
  # draw and save top down & up gene
  write.csv(data_down_gene_top30,file.path(output_dir,"deg_top30_down_gene_pvalueue.csv"))
  write.csv(data_up_gene_top30,file.path(output_dir,"deg_top30_up_gene_pvalueue.csv"))
  # draw down gene
  data_up_gene_top30_fc <- head(data_up_gene[order(-data_up_gene$log2FoldChange, -data_up_gene$pvalue), ],30)
  write.csv(data_up_gene_top30_fc,file.path(output_dir,"deg_top30_up_gene_fc.csv"))
  data_down_gene_top30_fc <- head(data_down_gene[order(data_down_gene$log2FoldChange, data_down_gene$pvalue), ],30)
  write.csv(data_down_gene_top30_fc,file.path(output_dir,"deg_top30_down_gene_fc.csv"))
  # 
  p3 <- draw_volcano_up_down(data,data_up_gene_top30_fc,data_down_gene_top30_fc,line_data_1,line_data_2)
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(FC).pdf")),plot = p3,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(FC).png")),device = "png",plot = p3,width = 9,height = 10,units = "cm",dpi = 1000)
  # DRAW EVERY one protein
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of up fc gene DO!!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_up_gene_top30_fc$SYMBOL), clear = FALSE, width = 100
  )
  name <- "up_fc"
  for(i in data_up_gene_top30_fc$SYMBOL){
    pb$tick()
    draw_volcano_specail(deg_result,i,name,figure_special_up_dir_fc)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of up fc gene DONE!!!!")))
  cat("\n")
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of down fc gene DO!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_down_gene_top30_fc$SYMBOL), clear = FALSE, width = 100
  )
  name <-"down_fc" 
  for(i in data_down_gene_top30_fc$SYMBOL){
    pb$tick()
    draw_volcano_specail(deg_result,i,name,figure_special_down_dir_fc)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of down fc gene DONE!!!!")))
  cat("\n")
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of down P-value  gene DO!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_down_gene_top30$SYMBOL), clear = FALSE, width = 100
  )
  name <- "down_pvalueue"
  for(i in data_down_gene_top30$SYMBOL){
    pb$tick()
    draw_volcano_specail(deg_result,i,name,figure_special_down_dir_pvalueue)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of down P-value  gene DONE!!!")))
  cat("\n")
  # draw up gene
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of up P-value gene DO!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_up_gene_top30$SYMBOL), clear = FALSE, width = 100
  )
  name <- "up_pvalueue"
  for(i in data_up_gene_top30$SYMBOL){
    pb$tick()
    draw_volcano_specail(deg_result,i,name,figure_special_up_dir_pvalueue)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of up P-value gene DONE!!!")))
  cat("\n")
}
# Function22: draw volcano style 2
draw_volcano_type_2_padj <- function(deg_result_add_lable){
  data <- deg_result
  data$FC <- 2^(data$log2FoldChange)
  data$log10 <- -log10(data$padj)
  pvalue = 0.05
  log2FC = 1
  f <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_1 <- data.frame(x = inputx + log2FC, y = y)
    return(dff_1)
  }
  r <- function(x){
    inputx <- seq(0.02, x, by = 0.0001)
    y <- 1/(inputx) + (-log10(pvalue))
    dff_2 <- data.frame(x = -(inputx + log2FC), y = y)
    return(dff_2)
  }
  line_data_1 <- f(5)
  line_data_2 <- r(5)
  data$curve_y <- case_when(
    data$log2FoldChange > 0 ~ 1/(data$log2FoldChange-log2FC) + (-log10(pvalue)),
    data$log2FoldChange <= 0 ~ 1/(-data$log2FoldChange-log2FC) + (-log10(pvalue)))
  data$group2 <- case_when(
    data$log10 > data$curve_y & data$log2FoldChange >= log2FC ~ 'Up-regulated',
    data$log10 > data$curve_y & data$log2FoldChange <= -log2FC ~ 'Down-regulated',
    TRUE ~ 'Non-regulated')
  head(data)
  data_up_gene <- subset(data,data$group2=="Up-regulated")
  data_up_gene_top10 <- head(data_up_gene[order(data_up_gene$padj),],10)
  data_up_gene_top30 <- head(data_up_gene[order(data_up_gene$padj),],30)
  data_down_gene <- subset(data,data$group2=="Down-regulated")
  data_down_gene_top10 <- head(data_down_gene[order(data_down_gene$padj),],10)
  data_down_gene_top30 <- head(data_down_gene[order(data_down_gene$padj),],30)
  # 1.5 save Volcano plot
  p1 <- draw_volcano_up_down_padj(data,data_up_gene_top10,data_down_gene_top10,line_data_1,line_data_2)
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top10 Volcano Plot type_2.pdf")),plot = p1,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top10 Volcano Plot type_2.png")),device = "png",plot = p1,width = 9,height = 10,units = "cm",dpi = 1000)
  #
  p2 <- draw_volcano_up_down_padj(data,data_up_gene_top30,data_down_gene_top30,line_data_1,line_data_2)
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(p-value).pdf")),plot = p2,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(p-value).png")),device = "png",plot = p2,width = 9,height = 10,units = "cm",dpi = 1000)
  
  # draw and save top down & up gene
  write.csv(data_down_gene_top30,file.path(output_dir,"deg_top30_down_gene_pvalueue.csv"))
  write.csv(data_up_gene_top30,file.path(output_dir,"deg_top30_up_gene_pvalueue.csv"))
  # draw down gene
  data_up_gene_top30_fc <- head(data_up_gene[order(-data_up_gene$log2FoldChange, -data_up_gene$padj), ],30)
  write.csv(data_up_gene_top30_fc,file.path(output_dir,"deg_top30_up_gene_fc.csv"))
  data_down_gene_top30_fc <- head(data_down_gene[order(data_down_gene$log2FoldChange, data_down_gene$padj), ],30)
  write.csv(data_down_gene_top30_fc,file.path(output_dir,"deg_top30_down_gene_fc.csv"))
  # 
  p3 <- draw_volcano_up_down_padj(data,data_up_gene_top30_fc,data_down_gene_top30_fc,line_data_1,line_data_2)
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(FC).pdf")),plot = p3,width = 9,height = 10,units = "cm")
  ggsave(file.path(figure_dir,paste0(get_parent_name()," top30  Volcano Plot type_2(FC).png")),device = "png",plot = p3,width = 9,height = 10,units = "cm",dpi = 1000)
  # DRAW EVERY one protein
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of up fc gene DO!!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_up_gene_top30_fc$SYMBOL), clear = FALSE, width = 100
  )
  name <- "up_fc"
  for(i in data_up_gene_top30_fc$SYMBOL){
    pb$tick()
    draw_volcano_specail_padj(deg_result,i,name,figure_special_up_dir_fc)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of up fc gene DONE!!!!")))
  cat("\n")
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of down fc gene DO!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_down_gene_top30_fc$SYMBOL), clear = FALSE, width = 100
  )
  name <- "down_fc"
  for(i in data_down_gene_top30_fc$SYMBOL){
    pb$tick()
    draw_volcano_specail_padj(deg_result,i,name,figure_special_down_dir_fc)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of down fc gene DONE!!!!")))
  cat("\n")
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of down Padj  gene DO!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_down_gene_top30$SYMBOL), clear = FALSE, width = 100
  )
  name <- "down_Padj"
  for(i in data_down_gene_top30$SYMBOL){
    pb$tick()
    draw_volcano_specail_padj(deg_result,i,name,figure_special_down_dir_pvalueue)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of down Padj  gene DONE!!!")))
  cat("\n")
  # draw up gene
  # process
  cat("\n")
  cat(bold(green("draw volcano plot of up Padj gene DO!!!")))
  cat("\n")
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total =length(data_up_gene_top30$SYMBOL), clear = FALSE, width = 100
  )
  name <- "up_Padj"
  for(i in data_up_gene_top30$SYMBOL){
    pb$tick()
    draw_volcano_specail_padj(deg_result,i,name,figure_special_up_dir_pvalueue)
  }
  cat("\n")
  cat(bold(green("draw volcano plot of up Padj gene DONE!!!")))
  cat("\n")
}
# Function23: draw volcano style 3
draw_volcano_type_3 <- function(deg_data,EXP_NAEE,figure_name,save_dir){
  # figure title
  if (length(strsplit(figure_name," ")[[1]])>6){
    prefix <- strsplit(figure_name,strsplit(figure_name," ")[[1]][4])[[1]][1]
    suffix <- strsplit(figure_name,strsplit(figure_name," ")[[1]][4])[[1]][2]
    new_title <- paste0(prefix,"\n",suffix)
    figure_name <- new_title
  }else{
    figure_name <- figure_name
  }
  #remove not suite condition for save file name
  if (grepl("\n",figure_name)){
    save_name <- gsub("\n","",figure_name)
  }else{
    save_name <- figure_name
  }
  # draw plot
  deg_result <- deg_data
  deg_result$log10 <- -log10(deg_result$pvalue)
  # ADD UP&DOWN&NO GENE TAG
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$pvalue < 0.05) & (deg_result$log2FoldChange > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$pvalue  < 0.05) & (deg_result$log2FoldChange < -1))] = "Down-regulated"
  # SORT BY PVALUE VALUE
  deg_result <- deg_result[order(deg_result$pvalue),]
  # GET NOT&UP&DOWN DATA
  non_deg_result <- subset(deg_result,deg_result$Group =="Non-significan")
  up_deg_result <- subset(deg_result,deg_result$Group =="Up-regulated")
  down_deg_result <- subset(deg_result,deg_result$Group =="Down-regulated")
  # GET TOP 15 pvalue GENE UP&DOWN
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),15)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),15)
  # get y_aes
  y_aes <- deg_result$log10
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  p <- ggplot(deg_result, aes(x = log2FoldChange, y = log10)) + 
    geom_point(data=non_deg_result,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,color="#C7C7C7",alpha=0.25) +
    geom_point(data=deg_result_up,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,fill="#e41749",alpha=0.5) +
    geom_point(data=up_deg_result,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,color="#e41749",alpha=0.4) +
    geom_point(data=deg_result_down,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,fill="#41b6e6",alpha=0.5) +
    geom_point(data=down_deg_result,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,color="#41b6e6",alpha=0.4) +
    geom_vline(xintercept=0,lty=2,col="black",lwd=0.1) +  
    geom_hline(yintercept = 1.3,lty=2,col="black",lwd=0.1) +  
    labs(x= bquote("RNA-seq " * log[2] * " fold change " * .(EXP_NAEE) * ":WT"),y= expression(paste(-log[10], "P-value")),title =paste0(figure_name," Volcano Plot")) + 
    geom_text_repel(data = deg_result_up,aes(log2FoldChange, log10, label= SYMBOL),size=1,colour="black",fontface="bold.italic",
                    segment.alpha = 0.5,segment.size = 0.15,segment.color = "black",min.segment.length=0,
                    box.padding=unit(0.2, "lines"),point.padding=unit(0, "lines"),force = 20,max.iter = 3e3,
                    xlim=c(0, 6),max.overlaps = 25,arrow=arrow(length = unit(0.02, "inches"))) + 
    geom_text_repel(data = deg_result_down,aes(log2FoldChange, log10, label= SYMBOL),size=1,colour="black",fontface="bold.italic",
                    segment.alpha =0.5,segment.size = 0.15,segment.color = "black",min.segment.length=0,
                    box.padding=unit(0.2, "lines"),point.padding=unit(0, "lines"),force = 20,max.iter = 3e3,
                    xlim=c(-6, 0),max.overlaps = 25,arrow=arrow(length = unit(0.02, "inches"))) + 
    annotate("text", x=-3.5, y=y_aes_value*0.97,size = 1.5,colour="#41b6e6", label= paste0("Down in ",EXP_NAEE))+
    geom_segment(aes(x = -2, y = y_aes_value*0.9, xend = -5.5, yend = y_aes_value*0.9),arrow = arrow(length = unit(0.15, "cm")),colour="#41b6e6",linewidth=0.4,
                 alpha=0.15,lineend="round",linejoin="round") +
    annotate("text", x=3.5, y=y_aes_value*0.97, size = 1.5,colour="#e41749",label= paste0("UP in ",EXP_NAEE) )+
    geom_segment(aes(x = 2, y = y_aes_value*0.9, xend = 5.5, yend = y_aes_value*0.9),arrow = arrow(length = unit(0.15, "cm")),colour="#e41749",linewidth=0.4,
                 alpha=0.15,lineend="round",linejoin="round") +
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),n.breaks = 8) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 10) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size =5,family="sans",face = "bold"),legend.position="none",
          legend.title=element_text(size =4,face = "bold",family="sans"), 
          text = element_text(size = 4,family="sans"),title = element_text(size = 4),           
          axis.title.x = element_text(size = 4,family="sans",face = "bold"),axis.title.y = element_text(size = 4,family="sans",face = "bold"),    
          axis.text.x = element_text(size = 4,family="sans"),axis.text.y = element_text(size = 4,family="sans"),
          panel.border = element_rect(color = "#606c70", fill = NA, size = 0.3),
          axis.line.x=element_line(linetype=1,color="#606c70",size=0.12),       
          axis.line.y=element_line(linetype=1,color="#606c70",size=0.12),
          axis.ticks.x=element_line(color="#606c70",size=0.12,lineend = 0.05),
          axis.ticks.length=unit(.08,"lines"),
          axis.ticks.y=element_line(color="#606c70",size=0.12,lineend = 0.05))    
  # 1.5 save Volcano plot
  ggsave(file.path(save_dir,paste0(save_name," Volcano Plot.pdf")),plot = p,width = 7,height = 5,units = "cm")
  ggsave(file.path(save_dir,paste0(save_name," Volcano Plot.png")),device = "png",plot = p,width = 7,height = 5,units = "cm",dpi = 1200)
}
# Function24: draw volcano style 3
draw_volcano_type_3_padj <- function(deg_data,EXP_NAEE,figure_name,save_dir){
  # figure title
  if (length(strsplit(figure_name," ")[[1]])>6){
    prefix <- strsplit(figure_name,strsplit(figure_name," ")[[1]][4])[[1]][1]
    suffix <- strsplit(figure_name,strsplit(figure_name," ")[[1]][4])[[1]][2]
    new_title <- paste0(prefix,"\n",suffix)
    figure_name <- new_title
  }else{
    figure_name <- figure_name
  }
  #remove not suite condition for save file name
  if (grepl("\n",figure_name)){
    save_name <- gsub("\n","",figure_name)
  }else{
    save_name <- figure_name
  }
  # draw plot
  deg_result <- deg_data
  deg_result$log10 <- -log10(deg_result$padj)
  # ADD UP&DOWN&NO GENE TAG
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$padj < 0.05) & (deg_result$log2FoldChange > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$padj < 0.05) & (deg_result$log2FoldChange < -1))] = "Down-regulated"
  # SORT BY PVALUE VALUE
  deg_result <- deg_result[order(deg_result$padj),]
  # GET NOT&UP&DOWN DATA
  non_deg_result <- subset(deg_result,deg_result$Group =="Non-significan")
  up_deg_result <- subset(deg_result,deg_result$Group =="Up-regulated")
  down_deg_result <- subset(deg_result,deg_result$Group =="Down-regulated")
  # GET TOP 15 pvalue GENE UP&DOWN
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),15)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),15)
  # get y_aes
  y_aes <- deg_result$log10
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  p <- ggplot(deg_result, aes(x = log2FoldChange, y = log10)) + 
    geom_point(data=non_deg_result,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,color="#C7C7C7",alpha=0.25) +
    geom_point(data=deg_result_up,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,fill="#e41749",alpha=0.5) +
    geom_point(data=up_deg_result,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,color="#e41749",alpha=0.4) +
    geom_point(data=deg_result_down,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,fill="#41b6e6",alpha=0.5) +
    geom_point(data=down_deg_result,aes(x = log2FoldChange, y = log10),size=0.02,shape = 21,color="#41b6e6",alpha=0.4) +
    geom_vline(xintercept=0,lty=2,col="black",lwd=0.1) +  
    geom_hline(yintercept = 1.3,lty=2,col="black",lwd=0.1) +  
    labs(x= bquote("RNA-seq " * log[2] * " fold change " * .(EXP_NAEE) * ":WT"),y= expression(paste(-log[10], "padj")),title =paste0(figure_name," Volcano Plot")) + 
    geom_text_repel(data = deg_result_up,aes(log2FoldChange, log10, label= SYMBOL),size=1,colour="black",fontface="bold.italic",
                    segment.alpha = 0.5,segment.size = 0.15,segment.color = "black",min.segment.length=0,
                    box.padding=unit(0.2, "lines"),point.padding=unit(0, "lines"),force = 20,max.iter = 3e3,
                    xlim=c(0, 6),max.overlaps = 25,arrow=arrow(length = unit(0.02, "inches"))) + 
    geom_text_repel(data = deg_result_down,aes(log2FoldChange, log10, label= SYMBOL),size=1,colour="black",fontface="bold.italic",
                    segment.alpha =0.5,segment.size = 0.15,segment.color = "black",min.segment.length=0,
                    box.padding=unit(0.2, "lines"),point.padding=unit(0, "lines"),force = 20,max.iter = 3e3,
                    xlim=c(-6, 0),max.overlaps = 25,arrow=arrow(length = unit(0.02, "inches"))) + 
    annotate("text", x=-3.5, y=y_aes_value*0.97,size = 1.5,colour="#41b6e6", label= paste0("Down in ",EXP_NAEE))+
    geom_segment(aes(x = -2, y = y_aes_value*0.9, xend = -5.5, yend = y_aes_value*0.9),arrow = arrow(length = unit(0.15, "cm")),colour="#41b6e6",linewidth=0.4,
                 alpha=0.15,lineend="round",linejoin="round") +
    annotate("text", x=3.5, y=y_aes_value*0.97, size = 1.5,colour="#e41749",label= paste0("UP in ",EXP_NAEE) )+
    geom_segment(aes(x = 2, y = y_aes_value*0.9, xend = 5.5, yend = y_aes_value*0.9),arrow = arrow(length = unit(0.15, "cm")),colour="#e41749",linewidth=0.4,
                 alpha=0.15,lineend="round",linejoin="round") +
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),n.breaks = 8) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 10) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size =5,family="sans",face = "bold"),legend.position="none",
          legend.title=element_text(size =4,face = "bold",family="sans"), 
          text = element_text(size = 4,family="sans"),title = element_text(size = 4),           
          axis.title.x = element_text(size = 4,family="sans",face = "bold"),axis.title.y = element_text(size = 4,family="sans",face = "bold"),    
          axis.text.x = element_text(size = 4,family="sans"),axis.text.y = element_text(size = 4,family="sans"),
          panel.border = element_rect(color = "#606c70", fill = NA, size = 0.3),
          axis.line.x=element_line(linetype=1,color="#606c70",size=0.12),       
          axis.line.y=element_line(linetype=1,color="#606c70",size=0.12),
          axis.ticks.x=element_line(color="#606c70",size=0.12,lineend = 0.05),
          axis.ticks.length=unit(.08,"lines"),
          axis.ticks.y=element_line(color="#606c70",size=0.12,lineend = 0.05))    
  # 1.5 save Volcano plot
  ggsave(file.path(save_dir,paste0(save_name," Volcano Plot.pdf")),plot = p,width = 7,height = 5,units = "cm")
  ggsave(file.path(save_dir,paste0(save_name," Volcano Plot.png")),device = "png",plot = p,width = 7,height = 5,units = "cm",dpi = 1200)
}
# Function25: draw volcano
draw_Volcano <- function(deg_result){
  if(min(na.omit(deg_result$pvalue)) == 0){
    print_color_note_type2_warring("P-value have zero !!! draw Volcano at padj levels !!!")
    print_color_note_type1("draw Volcano at padj levels DO!!!!!")
    # draw all result
    deg_result_volcano <- draw_volcano_padj(deg_result)
    deg_result_volcano_type_2 <- draw_volcano_type_2_padj(deg_result)
    # draw special gene volcano
    draw_volcano_specail_padj(deg_result,gene_1,"special",figure_special_dir)
    # draw Volcano plot paper
    draw_volcano_type_3_padj(deg_result,EXP_NAEE,figure_name,figure_dir)
    print_color_note_type3("draw Volcano at padj levels DONE!!!!!")
  }else{
    print_color_note_type1("draw Volcano at P-VALUE levels DO!!!!!")
    # draw all result
    deg_result_volcano <- draw_volcano(deg_result)
    deg_result_volcano_type_2 <- draw_volcano_type_2(deg_result)
    # draw special gene volcano
    draw_volcano_specail(deg_result,gene_1,"special",figure_special_dir)
    # draw Volcano plot paper
    draw_volcano_type_3(deg_result,EXP_NAEE,figure_name,figure_dir)
    print_color_note_type3("draw Volcano at P-VALUE levels DONE!!!!!")
  }
}
# Function26: draw activate volcano
draw_volcano_activate <- function(deg,i){
  # deg <- deg_result
  # draw plot
  deg_result <- deg
  deg_result$log10 <- -log10(deg_result$pvalue)
  # ADD UP&DOWN&NO GENE TAG
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$pvalue < 0.05) & (deg_result$log2FoldChange > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$pvalue  < 0.05) & (deg_result$log2FoldChange < -1))] = "Down-regulated"
  # SORT BY pvalue VALUE
  deg_result <- deg_result[order(deg_result$pvalue),]
  # GET NOT&UP&DOWN DATA
  non_deg_result <- subset(deg_result,deg_result$Group =="Non-significan")
  up_deg_result <- subset(deg_result,deg_result$Group =="Up-regulated")
  down_deg_result <- subset(deg_result,deg_result$Group =="Down-regulated")
  # GET TOP 15 pvalue GENE UP&DOWN
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),15)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),15)
  # get y_aes
  y_aes <- deg_result$log10
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  pal <- c("#F67280","#2daf94","gray")
  deg_result$Group <- factor(deg_result$Group,level=c("Up-regulated","Down-regulated","Non-significan"))
  plot <- plot_ly(deg_result,
                  x = ~log2FoldChange,
                  y = ~log10,
                  text = ~SYMBOL, 
                  color = ~Group,
                  colors = pal,
                  type = 'scatter',
                  mode = 'markers',
                  alpha = 0.5)
  # add Volcano Plot part
  plot <- layout(plot, 
                 xaxis = list(range = c(-1.2 * max(deg_result$log2FoldChange), 1.2 * max(deg_result$log2FoldChange)),title = "FoldChange(log2)"),
                 yaxis = list(range = c(0, max(-log10(deg_result$pvalue))),title = "-log10(P-value)"),
                 showlegend = TRUE,
                 legend = list(x = 1.05, y= 0))
  saveWidget(plot, file = file.path(figure_volcano_activate_dir,paste0(i," Volcano plot.html")), selfcontained = TRUE)
}
# Function27: draw activate volcano shiny
deg_result_shiny <- function(){
  # deg <- deg_result
  # draw plot
  deg_result$log10 <- -log10(deg_result$pvalue)
  # ADD UP&DOWN&NO GENE TAG
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$pvalue < 0.05) & (deg_result$log2FoldChange > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$pvalue  < 0.05) & (deg_result$log2FoldChange < -1))] = "Down-regulated"
  # SORT BY pvalue VALUE
  deg_result <- deg_result[order(deg_result$pvalue),]
  # GET NOT&UP&DOWN DATA
  non_deg_result <- subset(deg_result,deg_result$Group =="Non-significan")
  up_deg_result <- subset(deg_result,deg_result$Group =="Up-regulated")
  down_deg_result <- subset(deg_result,deg_result$Group =="Down-regulated")
  # GET TOP 15 pvalue GENE UP&DOWN
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),15)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),15)
  # get y_aes
  y_aes <- deg_result$log10
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$log2FoldChange
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  pal <- c("#F67280","#2daf94","gray")
  deg_result$Group <- factor(deg_result$Group,level=c("Up-regulated","Down-regulated","Non-significan"))
  plot <- plot_ly(deg_result,
                  x = ~log2FoldChange,
                  y = ~log10,
                  text = ~SYMBOL, 
                  color = ~Group,
                  colors = pal,
                  type = 'scatter',
                  mode = 'markers',
                  alpha = 0.5)
  # add Volcano Plot part
  plot <- layout(plot, 
                 xaxis = list(range = c(-1.2 * max(deg_result$log2FoldChange), 1.2 * max(deg_result$log2FoldChange)),title = "FoldChange(log2)"),
                 yaxis = list(range = c(0, max(-log10(deg_result$pvalue))),title = "-log10(P-value)"),
                 showlegend = TRUE,
                 legend = list(x = 1.05, y= 0))
  # use shiny build interactive Volcano Plot
  ui <- fluidPage(
    titlePanel("Dynamic Volcano Plot"),
    sidebarLayout(
      sidebarPanel(
        textInput("searchGene", "Search Gene", placeholder = "Enter gene name")
      ),
      mainPanel(
        plotlyOutput("volcanoPlot")
      )
    )
  )
  server <- function(input, output) {
    filteredPlot <- reactive({
      if (!is.null(input$searchGene) && input$searchGene != "") {
        selectedGene <- input$searchGene
        filteredData <- subset(deg_result, SYMBOL == selectedGene)
        if (nrow(filteredData) > 0) {
          plot %>%
            add_markers(data = filteredData,
                        x = ~log2FoldChange,
                        y = ~-log10(pvalue),
                        text = ~SYMBOL,
                        marker = list(color = 'red', size = 10))
        } else {
          plot
        }
      } else {
        plot
      }
    })
    output$volcanoPlot <- renderPlotly({
      filteredPlot()
    })
  }
  # run shiny
  shinyApp(ui = ui, server = server) 
}
# print run condition
print_color_note_type3("loading function done DONE!!!")
################################# PATH-5 #################################
# PATH-5  Set & Create run dir 
# print run condition
print_color_note_type1("Set & Create run dir DO!!!")
# set dir path
result_dir <- file.path(root_dir,output_name)
figure_dir <- file.path(result_dir,"figure")
output_dir <- file.path(result_dir,"output")
figure_special_dir <- file.path(figure_dir,"special")
figure_special_up_dir <- file.path(figure_dir,"up")
figure_special_up_dir_pvalueue <- file.path(figure_special_up_dir,"pvalueue")
figure_special_up_dir_fc <- file.path(figure_special_up_dir,"fc")
figure_special_down_dir <- file.path(figure_dir,"down")
figure_special_down_dir_pvalueue <- file.path(figure_special_down_dir,"pvalueue")
figure_special_down_dir_fc <- file.path(figure_special_down_dir,"fc")
figure_volcano_activate_dir <- file.path(figure_dir,"volcano_activate")
setwd(root_dir)
# dir list
list_dir <- c(result_dir,figure_dir,output_dir,figure_special_dir,figure_special_up_dir,figure_special_up_dir_pvalueue,
              figure_special_up_dir_fc,figure_special_down_dir,figure_special_down_dir_pvalueue,figure_special_down_dir_fc,
              figure_volcano_activate_dir)
# create folder
create_dir(list_dir)
# print run condition
print_color_note_type3("Set & Create run dir DONE!!!")
################################# PATH-6 #################################
# PATH-6  input count and sample infor file
# print run condition
print_color_note_type1("input count and sample infor file DO!!!")
# read count file 
sample_infor <- read.csv(sample_file_name_dir)
sample_name <- paste0(sample_infor$group_name,paste0("_",count_level,".read.count.tsv"))
# read frist file
df <- merge_count(count_dir,sample_name)
# delete gene_id row and set new row name
rownames(df) <- df$gene_id;df <- df[-1,];df <- df[,-1]
df_1 <- as.data.frame(lapply(df,as.numeric))
rownames(df_1) <- rownames(df)
colnames(df_1) <- sample_infor$group_name
rownames(sample_infor) <- sample_infor$group_name
# read sample infor file
sample_infor$condition <- factor(sample_infor$condition);rownames(sample_infor) <- sample_infor$group_name
# print run condition
print_color_note_type3("input count and sample infor file DONE!!!")
################################# PATH-7 #################################
# PATH-7  DEG
# print run condition
print_color_note_type1("DEG DO!!!")
# check count and infor Consistency
all(rownames(sample_infor) %in% colnames(df_1))
all(rownames(sample_infor) == colnames(df_1))
# print deg condition
print_color_note_type1("DEG analysis by DEseq2 DO !!!")
# deg
deg_result <- DEG(df_1,sample_infor,ref_name,paste0(get_parent_name(),"_DEG_output_result.csv"))
# print deg condition
print_color_note_type3("DEG analysis by DEseq2 DONE !!!")
# EXTERT gene list report
get_up_gene_list_padj <- get_up_gene_padj(deg_result)
get_down_gene_list_padj <- get_down_gene_padj(deg_result)
get_up_gene_list_pvalueue <- get_up_gene_pvalueue(deg_result)
get_down_gene_list_pvalueue <- get_down_gene_pvalueue(deg_result)
# get DEG result
pvalueue_up <- dim(get_up_gene_list_pvalueue)[1]
padj_up <- dim(get_up_gene_list_padj)[1]
pvalueue_down <- dim(get_down_gene_list_pvalueue)[1]
padj_down <- dim(get_down_gene_list_padj)[1]
# DEG infor
DEG_infor <- c(pvalueue_up,pvalueue_down,padj_up,padj_down)
DEG_infor <- as.data.frame(DEG_infor)
rownames(DEG_infor) <- c("pvalueue_up","pvalueue_down","padj_up","padj_down")
# T CONVERT  DATA
DEG_infor <- t(DEG_infor)
write.csv(DEG_infor,file.path(output_dir,"DEG_result.csv"))
# print run condition
print_color_note_type3("DEG DONE!!!")
################################# PATH-8 #################################
# PATH-8  PCA analysis
# print run condition
print_color_note_type1("PCA analysis DO!!!")
# PCA analysis
plot_PCA(df_1,sample_infor,paste0(get_parent_name(),"_PCA_plot"))
# print run condition
print_color_note_type2("PCA analysis DONE!!!")
################################# PATH-9 #################################
# PATH-9  heatmap corr
# print run condition
print_color_note_type1("heatmap corr DO!!!")
# draw corr heatmap
plot_heatmap_corr(deg_result)
# print run condition
print_color_note_type3("heatmap corr DONE!!!")
################################ PATH-10 #################################
# PATH-10  Volcano Plot
# print run condition
print_color_note_type1("Volcano Plot DO!!!")
# draw normal volcano
draw_Volcano(deg_result)
# draw activate volcano
draw_volcano_activate(deg_result,EXP_NAEE)
# open shiny volcano
deg_result_shiny()
# print run condition
print_color_note_type3("Volcano Plot DONE!!!")
############################### PATH-END #################################
