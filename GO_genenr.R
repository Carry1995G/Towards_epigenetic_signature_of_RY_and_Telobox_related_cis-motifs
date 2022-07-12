##GO-Analysis Visualization for multiple organisms#

library(xlsx)
library(dplyr)
library(pgirmess)
library(tidyr)
library(stringr)

# set the working directory where the tables to use are located
setwd("~/Studium/Master/Projektmodul Turck/GO-Analyse")

Genenumbers<-c()

###A thaliana####
setwd("Z:/dep_coupland/grp_turck/carina/Athaliana/GO-Analysis")

K27posGenes<-scan(
  "Z:/dep_coupland/grp_turck/carina/Athaliana/Chrom_FeatureOverlap/H3K27me3+_genes.txt",
  quiet=TRUE,
  what=""
)

RY=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telobox=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telolike=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")

Gene <- list(K27posGenes)
TargetGeneClusters=list.files(pattern="\\GenesinClusters\\.txt$", recursive=TRUE)

#### Calculations ####
Background<-lengths(Gene)  
Organism <- c("Athaliana")

for (i in 1:length(TargetGeneClusters)) {
  
  unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster")) 

  Background<-lengths(Gene)
  splitfilename=strsplit(TargetGeneClusters[[i]], '[-_|.]')[[1]]
  teloclustername=paste0(splitfilename[1],"_",splitfilename[2])
  
  if  ( splitfilename[1] == "RY" ) {
    TargetGenes_RY <- nrow(unfilteredTargetGenes %>% filter(cluster %in% RY)) 
  } else if ( splitfilename[1]  ==  "Telobox" | splitfilename[1] == "telobox") {
    TargetGenes_TB <- nrow(unfilteredTargetGenes %>% filter(cluster %in% Telobox)) 
  } else if ( splitfilename[1]  ==  "Telolike" | splitfilename[1] == "telolike" ) {
    TargetGenes_TL <- nrow(unfilteredTargetGenes %>% filter(cluster %in% Telolike)) 
  }
}
newGennumbers<-data.frame(Organism,Background,TargetGenes_RY,TargetGenes_TB,TargetGenes_TL)

Genenumbers<- rbind(Genenumbers,newGennumbers)

###File####

setwd("~/Studium/Master/Master thesis/Projektmodul Turck/GO-Analyse")

write.xlsx(Genenumbers, file = "GO_Genenumbers.xlsx")
