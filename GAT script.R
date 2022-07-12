#install.packages("ggsignif")
#install.packages("Rtools")
#install.packages("ggsci")
#install.packages("ggpubr")
#install.packages("pheatmap")
#install.packages("pals")

library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(reshape2)
library(dplyr)
library("ggsci")
library(pheatmap)
library("ggsignif")
library(readr)
library(gridExtra)
library(tidyr)
library(pals)



#setwd("~/netscratch/dep_coupland/grp_turck/carina/Tomato/GAT")
setwd("~/Studium/Master/Master thesis/Gat-Analysis")


Motifs=c("RY","Telobox","Telolike")
#Files=c("Telobox_GAT_result_Homer-Annot.txt","Telolike_GAT_result_Homer-Annot.txt", "RY_GAT_result_Homer-Annot.txt")

####Mpolymorpha####
out <- c("CDS", "mRNA", "primary_transcript","transcript","pseudogene", "tRNA", "rRNA", "lnc_RNA", "miRNA", "snRNA", "snoRNA", "cDNA_match")

for (i in 1:length(Motifs)) {
    File <- paste(Motifs[i],"_Marchantia.bed_GAT_result.txt",sep = "") 
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3)) %>%
    mutate(H3K27me3=ifelse(H3K27me3=="H3K27me+","H3K27me3+","H3K27me3-"))
  GAT <- GAT %>% filter(annotation != out) %>%
    mutate(annotation=replace(annotation, annotation=="five_prime_UTR", "5'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="three_prime_UTR", "3'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="1kbpromoter", "1kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="3kbpromoter", "3kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="100bpromoter", "100bp promoter"))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kb promoter" , "1kb promoter", "100bp promoter","TSS", "5'UTR", "exon", "intron" ,"gene", "3'UTR", "TTS","downprox","intergenic"))
  
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs")) +
    theme(text=element_text(size = 14))+
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Mpolymorpha_GAT_", Motifs[i],".png"), plot = a)
}

####Gmax######
out <- c("CDS", "mRNA", "primary_transcript","transcript","pseudogene", "tRNA", "rRNA", "lnc_RNA", "miRNA", "snRNA", "snoRNA", "cDNA_match")

for (i in 1:length(Motifs)) {
  
  File <- paste("Gmax_",Motifs[i],"GAT_result.txt",sep = "") 
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kbpromoter" , "1kbpromoter", "100bpromoter","TSS", "five_prime_UTR", "exon", "intron" ,"gene", "TTS","downprox", "intergenic"))
  GAT <- GAT %>% filter(annotation != out)
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs in G. max")) +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Gmax_GAT_", Motifs[i],".png"), plot = a)
}




####Bdistachyon####
setwd("Z:/dep_coupland/grp_turck/carina/Bdistachyon/GAT")

out <- c("CDS", "mRNA", "primary_transcript","transcript","pseudogene", "tRNA", "rRNA", "lnc_RNA", "miRNA", "snRNA", "snoRNA", "cDNA_match")

for (i in 1:length(Motifs)) {
  
  #File <- paste(Motifs[i],"GAT_result_Bdist.txt",sep = "") 
  File <- paste(Motifs[i],"GAT_result.txt",sep = "") 
  
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3)) %>% 
    mutate(annotation=replace(annotation, annotation=="five_prime_UTR", "5'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="three_prime_UTR", "3'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="1kbpromoter", "1kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="3kbpromoter", "3kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="100bpromoter", "100bp promoter"))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kb promoter" , "1kb promoter", "100bp promoter","TSS", "5'UTR", "exon", "intron" ,"gene", "3'UTR", "TTS","downprox","intergenic"))
  GAT <- GAT %>% filter(annotation != out)
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs")) +
    theme(text=element_text(size = 14))+
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Bdistachyon_GAT_", Motifs[i],".png"), plot = a,path="Z:/dep_coupland/grp_turck/carina/Bdistachyon/Plots")
}


####Medicago truncatula ####

out <- c("CDS", "mRNA", "primary_transcript","transcript","pseudogene", "tRNA", "rRNA", "lnc_RNA", "miRNA", "snRNA", "snoRNA", "cDNA_match", "inverted repeat",
         "PR_tract", "pre_miRNA")

for (i in 1:length(Motifs)) {
  
  File <- paste("Mtruncatula_",Motifs[i],"GAT_result.txt",sep = "") 
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3))
  GAT <- GAT %>% 
    mutate(annotation=replace(annotation, annotation=="five_prime_UTR", "5'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="three_prime_UTR", "3'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="1kbpromoter", "1kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="3kbpromoter", "3kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="100bpromoter", "100bp promoter"))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kb promoter" , "1kb promoter", "100bp promoter","TSS", "5'UTR", "exon", "intron" ,"gene", "3'UTR", "TTS","downprox","intergenic"))
  GAT <- GAT %>% filter(annotation != out)
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    theme(text=element_text(size = 14))+
    ggtitle(paste0(Motifs[i]," motifs")) +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Mtruncatula_GAT_", Motifs[i],".png"), plot = a)
}

####Tomato####
setwd("Z:/dep_coupland/grp_turck/carina/Tomato/GAT/GAT_Fast")

#GAT without any changes to workspace (except for without contigs)
for (i in 1:length(Motifs)) {
  File <- paste(Motifs[i],"_GAT_minimum.txt",sep = "") 
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3)) %>%
    mutate(H3K27me3=ifelse(H3K27me3=="H3K27me+","H3K27me3+","H3K27me3-"))
  GAT <- GAT %>% #filter(annotation != out) %>%
    mutate(annotation=replace(annotation, annotation=="five_prime_UTR", "5'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="three_prime_UTR", "3'UTR")) %>%
    mutate(annotation=replace(annotation, annotation=="1kbpromoter", "1kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="3kbpromoter", "3kb promoter")) %>%
    mutate(annotation=replace(annotation, annotation=="100bpromoter", "100bp promoter"))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kb promoter" , "1kb promoter", "100bp promoter","TSS", "5'UTR", "exon", "intron" ,"gene", "3'UTR", "TTS","downprox","intergenic"))
  
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs")) +
    theme(text=element_text(size = 14))+
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Tomato_minGAT_", Motifs[i],".png"), plot = a, path="Z:/dep_coupland/grp_turck/carina/Tomato/Plots")
}

#####Tomato plots different Chromosomes#####
for (i in 1:length(Motifs)) {
  pdf(file = paste0("Slycopersicum_GAT_", Motifs[i],"_nosmallintergenic.pdf"), onefile=TRUE)
  
  for (nr in 1:12) {
    File <- paste(Motifs[i],"GAT_result_CHR", nr, "_nosmallintergenic.txt",sep = "")

    #Telolike_-GAT_result_CHROM13
    #Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)

    GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
      mutate(H3K27me3=paste0("H3K",H3K27me3)) %>% 
      mutate(H3K27me3=ifelse(H3K27me3=="H3K27me+","H3K27me3+","H3K27me3-"))
    GAT <- GAT %>% 
      mutate(annotation=replace(annotation, annotation=="five_prime_UTR", "5'UTR")) %>%
      mutate(annotation=replace(annotation, annotation=="three_prime_UTR", "3'UTR")) %>%
      mutate(annotation=replace(annotation, annotation=="1kbpromoter", "1kb promoter")) %>%
      mutate(annotation=replace(annotation, annotation=="3kbpromoter", "3kb promoter")) %>%
      mutate(annotation=replace(annotation, annotation=="100bpromoter", "100bp promoter"))
    GAT$annotation <- factor(GAT$annotation,levels = c("3kbpromoter" , "1kbpromoter", "100bpromoter","TSS", "five_prime_UTR", "exon", "intron" ,"gene", "three_prime_UTR", "TTS","downprox","intergenic"))
    
    a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
      geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
      scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
      ggtitle(paste0(Motifs[i],"motifs in S. lycopersicum Chr ",nr-1)) +
      geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
    
    print(a)
  }
  dev.off()
}

#####Tomato plots for 1st chrom vs last chroms ####

for (i in 1:length(Motifs)) {
  Filecontigs <- paste(Motifs[i],"_GAT_result_contigs.txt",sep = "") 
  FileCHROMs <- paste(Motifs[i],"_GAT_result_CHROMs.txt",sep = "") 
  
  GAT_resultcontigs <- read_delim(Filecontigs, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (Chrom ="contigs") %>% select(-1)
  GAT_resultCHROMs <- read_delim(FileCHROMs, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (Chrom ="chromosomes") %>% select(-1)
  GAT_result <- rbind(GAT_resultcontigs,GAT_resultCHROMs)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kbpromoter" , "1kbpromoter", "100bpromoter","TSS", "five_prime_UTR", "exon", "intron" ,"gene", "three_prime_UTR", "TTS","downprox","intergenic"))
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ Chrom + H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs in S. lycopersicum")) +
    geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Slycopersicum_contigs+Chrs_GAT_", Motifs[i],".png"), plot = a)
}

#####no repeats########
for (i in 1:length(Motifs)) {

  File <- paste(Motifs[i],"GAT_result_repeats+BL4.txt",sep = "") 

  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)

  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kbpromoter" , "1kbpromoter", "100bpromoter","TSS", "five_prime_UTR", "exon", "intron" ,"gene", "three_prime_UTR", "TTS","downprox","intergenic"))
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs in S. lycopersicum")) +
    geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Slycopersicum_repeat+BL4_GAT_", Motifs[i],".png"), plot = a)
}

#####no intergenic##########

Regions <-c("nointergenic", "notelorpts+rpts+intergenic", "norpts+intergenic")
for (r in Regions) {
for (i in 1:length(Motifs)) {
  
  File <- paste(Motifs[i],"GAT_result_",r,".txt",sep = "") 
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3)) %>% filter(annotation!='intergenic')
  GAT$annotation <- factor(GAT$annotation,levels = c("3kbpromoter" , "1kbpromoter", "100bpromoter","TSS", "five_prime_UTR", "exon", "intron" ,"gene", "three_prime_UTR", "TTS","downprox","intergenic"))
  
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs in S. lycopersicum")) +
    geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Slycopersicum_",r,"_", Motifs[i],".png"), plot = a)
}
}


#####no intergenic, but less####
for (i in 1:length(Motifs)) {
  
  File <- paste(Motifs[i],"GAT_result_nosmallintergenic.txt",sep = "") 
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)
  
  GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
    mutate(H3K27me3=paste0("H3K",H3K27me3)) %>% filter(annotation!='intergenic') %>% mutate(H3K27me3=ifelse(H3K27me3=="H3K27me+","H3K27me3+","H3K27me3-"))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kbpromoter" , "1kbpromoter", "100bpromoter","TSS", "five_prime_UTR", "exon", "intron" ,"gene", "three_prime_UTR", "TTS","downprox","intergenic"))
  
  
  a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
    geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
    ggtitle(paste0(Motifs[i]," motifs in S. lycopersicum")) +
    geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  ggsave (paste0("Slycopersicum_lessintergenic", Motifs[i],".png"), plot = a)
}
#####Tomato plots for different chromosomes + - strain####
for (i in 1:length(Motifs)) {
  pdf(file = paste0("GAT_", Motifs[i],".pdf"), onefile=TRUE)
  
  for (nr in 1:13) {
    Fileneg <- paste(Motifs[i],"_-GAT_result_CHROM", nr, ".txt",sep = "")
    Filepos <- paste(Motifs[i],"_+GAT_result_CHROM", nr,".txt",sep = "")
    
    #Telolike_-GAT_result_CHROM13
    #Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    GAT_resultneg <- read_delim(Fileneg, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (strain ="- strain") %>% select(-1)
    GAT_resultpos <- read_delim(Filepos, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (strain ="+ strain") %>% select(-1)
    GAT_result <- rbind(GAT_resultneg,GAT_resultpos)
    
    GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
      mutate(H3K27me3=paste0("H3K",H3K27me3))
    
    a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
      geom_col() + facet_wrap(~ strain + H3K27me3) + scale_fill_brewer(palette = "Paired") +
      scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
      ggtitle(paste0("CHR",nr)) +
      geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
    
    print(a)
  }
  dev.off()
}



#####Tomato test ####
for (i in 1:length(Motifs)) {
  pdf(file = paste0("GAT_", Motifs[i],".pdf"), onefile=TRUE)
  
  for (nr in 1:13) {
Fileneg <- paste(Motifs[i],"_-GAT_result_CHROM", nr, ".txt",sep = "")
Filepos <- paste(Motifs[i],"_+GAT_result_CHROM", nr,".txt",sep = "")

#Telolike_-GAT_result_CHROM13
  #Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
  GAT_resultneg <- read_delim(Fileneg, "\t", escape_double = FALSE, trim_ws = TRUE)
  GAT_resultpos <- read_delim(Filepos, "\t", escape_double = FALSE, trim_ws = TRUE)
  

  GAT <- GAT_resultneg %>% dplyr::filter(annotation %in% c("3kbpromoter_K27+" , "1kbpromoter_K27+", "100bpromoter_K27+", "proxp_K27+","TSS_K27+", "five_prime_UTR_K27+", "exon_K27+","CDS_K27+", "intron_K27+","three_prime_UTR_K27+", "TTS_K27+","downprox_K27+","intergenic_K27+"))
  GAT$annotation <- factor(GAT$annotation,levels = c("3kbpromoter_K27+" , "1kbpromoter_K27+", "100bpromoter_K27+", "proxp_K27+","TSS_K27+", "five_prime_UTR_K27+", "exon_K27+","CDS_K27+", "intron_K27+","three_prime_UTR_K27+", "TTS_K27+","downprox_K27+","intergenic_K27+"))
  a<-ggplot(GAT)+geom_col(aes(annotation,l2fold,fill=annotation),  position = "dodge") + 
    theme_bw(base_size=12) + scale_fill_npg() + ggtitle(paste0("CHR",nr,"+ strain H3K27me3+ ")) +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  GAT1 <- GAT_resultneg %>% dplyr::filter(annotation %in% c("3kbpromoter_K27-" , "1kbpromoter_K27-", "100bpromoter_K27-","proxp_K27-","TSS_K27-", "five_prime_UTR_K27-", "exon_K27-","CDS_K27-", "intron_K27-","three_prime_UTR_K27-", "TTS_K27-","downprox_K27-","intergenic_K27-"))
  GAT1$annotation <- factor(GAT1$annotation,levels = c("3kbpromoter_K27-" , "1kbpromoter_K27-", "100bpromoter_K27-","proxp_K27-","TSS_K27-", "five_prime_UTR_K27-", "exon_K27-","CDS_K27-", "intron_K27-","three_prime_UTR_K27-", "TTS_K27-","downprox_K27-","intergenic_K27-"))
  a2<-ggplot(GAT1)+geom_col(aes(annotation,l2fold), position = "dodge") + 
    theme_bw(base_size=12) + scale_fill_npg() + ggtitle(paste0("CHR",nr,"+ strain H3K27me3- ")) +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  
  GAT2 <- GAT_resultpos %>% dplyr::filter(annotation %in% c("3kbpromoter_K27+" , "1kbpromoter_K27+", "100bpromoter_K27+","proxp_K27+","TSS_K27+", "five_prime_UTR_K27+", "exon_K27+","CDS_K27+", "intron_K27+","three_prime_UTR_K27+", "TTS_K27+","downprox_K27+","intergenic_K27+"))
  GAT2$annotation <- factor(GAT2$annotation,levels = c("3kbpromoter_K27+" , "1kbpromoter_K27+", "100bpromoter_K27+","proxp_K27+","TSS_K27+", "five_prime_UTR_K27+", "exon_K27+","CDS_K27+", "intron_K27+","three_prime_UTR_K27+", "TTS_K27+","downprox_K27+","intergenic_K27+"))
  b<-ggplot(GAT2)+geom_col(aes(annotation,l2fold), position = "dodge") + 
    theme_bw(base_size=12) + scale_fill_npg() + ggtitle(paste0("CHR",nr,"- strain H3K27me3+ ")) +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  GAT3 <- GAT_resultpos %>% dplyr::filter(annotation %in% c("3kbpromoter_K27-" , "1kbpromoter_K27-", "100bpromoter_K27-","proxp_K27-","TSS_K27-", "five_prime_UTR_K27-", "exon_K27-","CDS_K27-", "intron_K27-","three_prime_UTR_K27-", "TTS_K27-","downprox_K27-","intergenic_K27-"))
  GAT3$annotation <- factor(GAT3$annotation,levels = c("3kbpromoter_K27-" , "1kbpromoter_K27-", "100bpromoter_K27-","proxp_K27-","TSS_K27-", "five_prime_UTR_K27-", "exon_K27-","CDS_K27-", "intron_K27-","three_prime_UTR_K27-", "TTS_K27-","downprox_K27-","intergenic_K27-"))
  b2<-ggplot(GAT3)+geom_col(aes(annotation,l2fold), position = "dodge") + 
    theme_bw(base_size=12) + scale_fill_npg() + ggtitle(paste0("CHR",nr," - strain H3K27me3- ")) + 
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  grid.arrange(arrangeGrob(a,a2,ncol=2), arrangeGrob(b,b2,ncol=2))
  }
  dev.off()
}

######### Tomato - + strain 

pdf(file = "GAT_Slycopersicum_+-strain.pdf", onefile=TRUE)
for (i in 1:length(Motifs)) {
  Fileneg <- paste(Motifs[i],"_-GAT_result_GFF-Annot.txt",sep = "")
  Filepos <- paste(Motifs[i],"_+GAT_result_GFF-Annot.txt",sep = "")

#Telolike_-GAT_result_CHROM13
#Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GAT_resultneg <- read_delim(Fileneg, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (strain ="- strain") %>% select(-1)
GAT_resultpos <- read_delim(Filepos, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (strain ="+ strain") %>% select(-1)
GAT_result <- rbind(GAT_resultneg,GAT_resultpos)

GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
  mutate(H3K27me3=paste0("H3K",H3K27me3))

a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
  geom_col() + facet_wrap(~ strain + H3K27me3) + scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
  ggtitle(paste0(Motifs[i]," Motif in S. lycopersicum")) +
  geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
print(a)
}
dev.off()
  

####Athaliana####
#####Athaliana Plots + - strain ####

Fileneg <- "RY_-_GAT_result.txt"
Filepos <- "RY_+_GAT_result.txt"

#Telolike_-GAT_result_CHROM13
#Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GAT_resultneg <- read_delim(Fileneg, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (strain ="- strain") %>% select(-1)
GAT_resultpos <- read_delim(Filepos, "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate (strain ="+ strain") %>% select(-1)
GAT_result <- rbind(GAT_resultneg,GAT_resultpos)

GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
  mutate(H3K27me3=paste0("H3K",H3K27me3))

a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
  geom_col() + facet_wrap(~ strain + H3K27me3) + scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
  ggtitle(paste0("RY Motif in A. thaliana ")) +
  geom_text(aes(annotation,l2fold,label = ifelse(pvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
pdf(file = "RY_Athaliana_GAT.pdf")
print(a)
dev.off()


###Athaliana Plot one strain ####

for (i in 1:length(Motifs)) {
  
  File <- paste("Athaliana_",Motifs[i],"_GAT_result.txt",sep = "") 
  
  GAT_result <- read_delim(File, "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(-1)

GAT <- GAT_result %>% separate(annotation, c("annotation", "H3K27me3"), "_K") %>%
  mutate(H3K27me3=paste0("H3K",H3K27me3)) %>% 
  mutate(H3K27me3=ifelse(H3K27me3=="H3K27+","H3K27me3+","H3K27me3-"))%>%
  mutate(annotation=replace(annotation, annotation=="1kbpro", "1kb promoter")) %>%
  mutate(annotation=replace(annotation, annotation=="five_prime_UTR", "5'UTR")) %>%
  mutate(annotation=replace(annotation, annotation=="three_prime_UTR", "3'UTR")) %>%
  mutate(annotation=replace(annotation, annotation=="3kbpro", "3kb promoter")) %>%
  mutate(annotation=replace(annotation, annotation=="proxp", "100bp promoter"))
GAT$annotation <- factor(GAT$annotation,levels = c("3kb promoter" , "1kb promoter", "100bp promoter","TSS", "5'UTR", "exon", "intron" ,"gene", "3'UTR", "TTS","downprox","intergenic"))


a<-ggplot(GAT, aes(x=annotation, y=l2fold, fill=annotation)) + 
  geom_col() + facet_wrap(~ H3K27me3) + scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") + 
  ggtitle(paste0(Motifs[i]," motifs")) +
  theme(text=element_text(size = 14))+
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
ggsave (paste0("Athaliana_GAT_", Motifs[i],".png"), plot = a)
}


###heatmap####

heat<-dcast(GAT, annotation ~ annotation)
rownames(heat) <- heat[,1]
heat[,1] <- NULL
pdf(file = "position_HM.pdf", width = 6, height = 4)
pheatmap(heat, cluster_cols=FALSE)
dev.off()



############## Plotting z-value ############
  az<-ggplot(GAT)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  bz<-ggplot(GAT2)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) + scale_fill_npg() +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  cz<-ggplot(GAT3)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) + scale_fill_npg() +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  pdf(file = paste0("GAT_z_0.05", Motifs[3],".pdf"), width = 12, height = 12)
  ggarrange(bz,cz, ncol = 1, nrow = 2)
  dev.off()





#############


quit(save="no")








Motifs=c("Telobox", "Telolike")
Files=c("Telobox_GAT_result.txt","Telolike_GAT_result.txt" )

for (i in 1:length(Motifs)) {
  
#Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GAT_result <- read_delim(Files[3], "\t", escape_double = FALSE, trim_ws = TRUE)


GAT <- GAT_result %>% dplyr::filter(annotation %in% c("TTS","intron" , "1kb-promoter", "exon", "TSS")) %>%
  mutate (Strain = )
GAT$annotation <- factor(GAT$annotation ,levels = c("TTS","intron" , "1kb-promoter", "exon", "TSS"))
a<-ggplot(GAT)+geom_col(aes(annotation,l2fold), position = "dodge") + theme_bw(base_size=12) +scale_fill_npg() +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.005, "*", "NS")),position = position_stack(vjust = 0.5))

GAT2 <- GAT_result %>% dplyr::filter(annotation %in% c("TTS_K27+","intron_K27+" , "1kb-promoter_K27+", "exon_K27+", "TSS_K27+"))
GAT2$annotation <- factor(GAT2$annotation ,levels = c("TTS_K27+","intron_K27+" , "1kb-promoter_K27+", "exon_K27+", "TSS_K27+"))
b<-ggplot(GAT2)+geom_col(aes(annotation,l2fold), position = "dodge") + theme_bw(base_size=12) +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.005, "*", "NS")),position = position_stack(vjust = 0.5))

GAT3 <- GAT_result %>% dplyr::filter(annotation %in% c("TTS_K27-","intron_K27-", "1kb-promoter_K27-", "exon_K27-", "TSS_K27-"))
GAT3$annotation <- factor(GAT3$annotation ,levels = c("TTS_K27-","intron_K27-", "1kb-promoter_K27-", "exon_K27-", "TSS_K27-"))
c<-ggplot(GAT3)+geom_col(aes(annotation,l2fold), position = "dodge") + theme_bw(base_size=12) +scale_fill_npg() +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.005, "*", "NS")),position = position_stack(vjust = 0.5))

pdf(file = paste0("GAT_", Motifs[3],".pdf"), width = 12, height = 12)
ggarrange(a,b,c, ncol = 1, nrow = 3)
dev.off()


az<-ggplot(GAT)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) +scale_fill_npg() +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.005, "*", "NS")),position = position_stack(vjust = 0.5))

bz<-ggplot(GAT2)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) +scale_fill_npg() +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.005, "*", "NS")),position = position_stack(vjust = 0.5))

cz<-ggplot(GAT3)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) +scale_fill_npg() +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.005, "*", "NS")),position = position_stack(vjust = 0.5))

pdf(file = paste0("GAT_z_", Motifs[3],".pdf"), width = 12, height = 12)
ggarrange(az,bz,cz, ncol = 1, nrow = 3)
dev.off()
}