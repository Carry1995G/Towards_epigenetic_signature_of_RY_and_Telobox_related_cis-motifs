#!/usr/bin/Rscript

#install.packages("RColorBrewer",repos= "http://cran.us.r-project.org")
#install.packages("ggplot2",repos= "http://cran.us.r-project.org")
#install.packages("MASS", repos= "http://cloud.r-project.org")
#install.packages("ggpubr", repos= "http://cloud.r-project.org", dependencies = TRUE)


#install.packages("ggsignif", repos= "http://cran.us.r-project.org")
#install.packages("Rtools", repos= "http://cran.us.r-project.org")
#install.packages("ggsci", repos= "http://cran.us.r-project.org")
#install.packages("pheatmap", repos= "http://cran.us.r-project.org")
#install.packages("readr", repos= "http://cran.us.r-project.org")

library(RColorBrewer)
library(ggplot2)

#library(ggsci)
#library(ggpubr)
#library(reshape2)
library(dplyr)
#library("ggsci")
#library(pheatmap)
#library("ggsignif")
#library("Rtools")
library(readr)

setwd("/netscratch/dep_coupland/grp_turck/carina/Tomato/GAT")


Motifs=c("Telobox", "Telolike", "RY")
Files=c("Telobox_GAT_result.txt","Telolike_GAT_result.txt", "RY_GAT_result.txt")

for (i in 1:length(Motifs)) {
  
  #Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
  GAT_result <- read_delim(Files[3], "\t", escape_double = FALSE, trim_ws = TRUE)
  
  
  GAT <- GAT_result %>% dplyr::filter(annotation %in% c("3kbpromoter_K27+" , "1kbpromoter_K27+", "100bpromoter_K27+", "proxp_K27+","TSS_K27+", "five_prime_UTR_K27+", "exon_K27+","CDS_K27+", "intron_K27+","three_prime_UTR_K27+", "TTS_K27+","downprox_K27+","intergenic_K27+"))
  GAT$annotation <- factor(GAT$annotation ,levels = c("TTS","intron" , "1kb-promoter", "3kb-promoter", "exon", "TSS"))
  a<-ggplot(GAT)+geom_col(aes(annotation,l2fold), position = "dodge") + theme_bw(base_size=12) +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  GAT2 <- GAT_result %>% dplyr::filter(annotation %in% c("3kbpromoter_K27+" , "1kbpromoter_K27+", "100bpromoter_K27+","proxp_K27+","TSS_K27+", "five_prime_UTR_K27+", "exon_K27+","CDS_K27+", "intron_K27+","three_prime_UTR_K27+", "TTS_K27+","downprox_K27+","intergenic_K27+"))
  GAT2$annotation <- factor(GAT2$annotation ,levels = c("3kbpromoter_K27+" , "1kbpromoter_K27+", "100bpromoter_K27+","proxp_K27+","TSS_K27+", "five_prime_UTR_K27+", "exon_K27+","CDS_K27+", "intron_K27+","three_prime_UTR_K27+", "TTS_K27+","downprox_K27+","intergenic_K27+"))
  b<-ggplot(GAT2)+geom_col(aes(annotation,l2fold), position = "dodge") + theme_bw(base_size=12) + scale_fill_npg() +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  GAT3 <- GAT_result %>% dplyr::filter(annotation %in% c("3kbpromoter_K27-" , "1kbpromoter_K27-", "100bpromoter_K27-","proxp_K27-","TSS_K27-", "five_prime_UTR_K27-", "exon_K27-","CDS_K27-", "intron_K27-","three_prime_UTR_K27-", "TTS_K27-","downprox_K27-","intergenic_K27-"))
  GAT3$annotation <- factor(GAT3$annotation ,levels = c("3kbpromoter_K27-" , "1kbpromoter_K27-", "100bpromoter_K27-","proxp_K27-","TSS_K27-", "five_prime_UTR_K27-", "exon_K27-","CDS_K27-", "intron_K27-","three_prime_UTR_K27-", "TTS_K27-","downprox_K27-","intergenic_K27-"))
  c<-ggplot(GAT3)+geom_col(aes(annotation,l2fold), position = "dodge") + theme_bw(base_size=12) + scale_fill_npg() +
    geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))
  
  #setwd("~/Studium/Master/Master thesis")
  pdf(file = paste0("GAT_0.05", Motifs[3],".pdf"), width = 12, height = 12)
  ggarrange(b,c, ncol = 1, nrow = 2)
  dev.off()
}

az<-ggplot(GAT)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))

bz<-ggplot(GAT2)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) + scale_fill_npg() +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))

cz<-ggplot(GAT3)+geom_col(aes(annotation,`-nan`), position = "dodge") + theme_bw(base_size=12) + scale_fill_npg() +
  geom_text(aes(annotation,l2fold,label = ifelse(qvalue < 0.05, "*", "NS")),position = position_stack(vjust = 0.5))

pdf(file = paste0("GAT_z_0.05", Motifs[3],".pdf"), width = 12, height = 12)
ggarrange(bz,cz, ncol = 1, nrow = 2)
dev.off()
}




#############


quit(save="no")


Motifs=c("Telobox", "Telolike")
Files=c("Telobox_GAT_result.txt","Telolike_GAT_result.txt" )

for (i in 1:length(Motifs)) {
  
  #Telolike_GAT_result <- read_delim("Telolike_GAT_result.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
  GAT_result <- read_delim(Files[3], "\t", escape_double = FALSE, trim_ws = TRUE)
  
  
  GAT <- GAT_result %>% dplyr::filter(annotation %in% c("TTS","intron" , "1kb-promoter", "exon", "TSS"))
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