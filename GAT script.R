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
library(ggsci)
library(pheatmap)
library(ggsignif)
library(readr)
library(gridExtra)
library(tidyr)
library(pals)

setwd("~")
Motifs=c("RY","Telobox","Telolike")

out <- c("CDS", "mRNA", "primary_transcript","transcript","pseudogene", "tRNA", "rRNA", "lnc_RNA", "miRNA", "snRNA", "snoRNA", "cDNA_match")

for (i in 1:length(Motifs)) {
    File <- paste(Motifs[i],"GAT_result.txt",sep = "") 
  
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
