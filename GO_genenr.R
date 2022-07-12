##GO-Analysis Visualization for multiple organisms#

# Requires the package 'ggplot2' (needs to be installed first)
# Load the ggplot2 package
library(ggplot2)
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
#### Marchantia polymorpha ####

setwd("Z:/dep_coupland/grp_turck/carina/Moss-Liverworts/Marchantia-chromatin/GO-Analyse")

geneID2GO <- (readMappings("GO-Annotation TopGO.map"))
names(geneID2GO)<-gsub("\\..*", "", names(geneID2GO)) 

#the Heatmap files contain wrong gene annotations: olyID=Mapoly0103s0044
#The correct annotations are included in the 10th Column of the Annotation files:
#chr1	531466	534441	olyID=Mapoly0103s0044	.	-	feature	gene	.	ID=Mp1g00430;MapolyID=Mapoly0103s0044;Name=Mp1g00430
#this code shrinks the annotation file to the needed columns (select) and separates the IDs into the needed one

translation <-read.table("Z:/dep_coupland/grp_turck/carina/Moss-Liverworts/Downloads/gene_correspondence_table_v51.txt",sep="\t",na.strings = "")
translationshort <- translation %>% drop_na %>% mutate (V1 = gsub("\\.[0-9]*$","",V1),
                                                        V2 = gsub("\\.[0-9]*$","",V2)) %>% distinct()
double <- translationshort[duplicated(translationshort$V2),] %>% distinct()

#not needed
#Genes4translate <- read.table("Z:/dep_coupland/grp_turck/carina/Moss-Liverworts/Annotation_Files/Moss_GFF3_genes.bed",sep="\t",header=FALSE) %>%
#dplyr::select(V4,V10) %>% filter(grepl('olyID=Mapoly', V4))  %>% 
#separate(V10, sep = "=|;", into = c("a", "Gene", "c")) %>% dplyr::select(V4,Gene) %>% mutate(V4 = gsub("olyID=","",V4))
#GenesandtranslatedGenes <- Genes4translate %>% left_join( translation, by=c("V4" = "V2"))
# mutate (newGenes = translation$V1[which(translation$V2==V4)])
#oldGO <- read.table("Z:/dep_coupland/grp_turck/carina/Moss-Liverworts/Downloads/gene.functions.txt",sep="\t", comment.char = "", header=TRUE)
#names(oldGO)[4]<-"description"
#rename(oldGO, Description = `description.if.any..GO.is.via.pfam2go.interpro`)
#oldGOtest <- oldGO %>% filter (IdType == "GO")
#oldGOtest <- oldGO %>% filter(grepl("GO",description))

Genes=scan(
  "genes.txt",
  quiet=TRUE,
  what=""
)

K27posGenes<-scan(
  "Z:/dep_coupland/grp_turck/carina/Moss-Liverworts/Marchantia-chromatin/FeatureOverlap/H3K27me3+_genes.txt",
  quiet=TRUE,
  what=""
)

RY=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telobox=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telolike=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")


#Telobox=c("cluster_3")
#Telolike=c("cluster_1")

Gene <- list(K27posGenes)
Genebackground<-c("K27+_Genes")
TargetGeneClusters=list.files(pattern="RY_Mpolymorpha_genesGOkmeans\\.txt_GenesinClusters\\.txt$", recursive=TRUE)
#TargetCDSClusters=list.files(pattern="\\CDSinClusters\\.txt$", recursive=TRUE)

#needed in loop:
#unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster"))  %>% rowwise() %>% #, col.names = c("gene", "cluster")) 
# mutate(gene=gsub(".*\\M","M",gene), 
#       gene=ifelse(gene %in% translationshort$V2[duplicated(translationshort$V2)],"DOUBLE",
#                  ifelse(gene %in% unique(translationshort$V2), translationshort$V1[which(translationshort$V2==gene)],
#                        gene)))


#### Calculations ####

Background<-lengths(Gene)  
Organism <- c("Mpolymorpha")

for (i in 1:length(TargetGeneClusters)) {
  
  unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster")) 
  
  
  splitfilename=strsplit(TargetGeneClusters[[i]], '[-_|.]')[[1]]
  teloclustername=paste0(splitfilename[1],"_",splitfilename[2])
  
  if  ( splitfilename[1] == "RY" ) {
    TargetGenes_RY <- nrow(unfilteredTargetGenes %>% filter(cluster %in% RY)) 
    Background_RY <-nrow(unfilteredTargetGenes)
    TargetGenes_TB <- 0 
    Background_TB <- 0 
    TargetGenes_TL <- 0 
    Background_TL <- 0 
  } else if ( splitfilename[1]  ==  "Telolike" | splitfilename[1] == "telolike" ) {
    TargetGenes_TL <- 0 
    Background_TL <- 0 
  }
}
newGennumbers<-data.frame(Organism,Background,TargetGenes_RY,TargetGenes_TB,TargetGenes_TL)
Genenumbers<- rbind(Genenumbers,newGennumbers)

###B distachyon####

setwd("Z:/dep_coupland/grp_turck/carina/Bdistachyon/GO-Analysis")


K27posGenes<-scan(
  "Z:/dep_coupland/grp_turck/carina/Bdistachyon/Chrom_FeatureOverlap/Genes/H3K27me3+genes.txt",
  quiet=TRUE,
  what=""
)

RY=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telobox=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telolike=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")

Gene <- list(K27posGenes)
Genebackground<-c("K27+_Genes")
TargetGeneClusters=list.files(pattern="\\genes_GenesinClusters\\.txt$", recursive=TRUE)
#### Calculations ####

Background<-lengths(Gene)
Organism <- c("Bdistachyon")

for (i in 1:length(TargetGeneClusters)) {
  
  unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster")) 
  
  
  splitfilename=strsplit(TargetGeneClusters[[i]], '[-_|.]')[[1]]
  teloclustername=paste0(splitfilename[1],"_",splitfilename[2])
  
  if  ( splitfilename[1] == "RY" ) {
    TargetGenes_RY <- nrow(unfilteredTargetGenes %>% filter(cluster %in% RY)) 
    Background_RY <-nrow(unfilteredTargetGenes)
  } else if ( splitfilename[1]  ==  "Telobox" | splitfilename[1] == "telobox") {
    TargetGenes_TB <- nrow(unfilteredTargetGenes %>% filter(cluster %in% Telobox)) 
    Background_TB <- nrow(unfilteredTargetGenes)
  } else if ( splitfilename[1]  ==  "Telolike" | splitfilename[1] == "telolike" ) {
    TargetGenes_TL <- nrow(unfilteredTargetGenes %>% filter(cluster %in% Telolike)) 
    Background_TL <-nrow(unfilteredTargetGenes)
  }
}
newGennumbers<-data.frame(Organism,Background,TargetGenes_RY,TargetGenes_TB,TargetGenes_TL)

Genenumbers<- rbind(Genenumbers,newGennumbers)
####Solanum lycopersicum####

setwd("Z:/dep_coupland/grp_turck/carina/Tomato/GO-Analysis")

#Gene files for TopGO: scan Gene txt files. Create vector with p-values and assign gene names to vector members 
#Genes=scan(
# "Z:/dep_coupland/grp_turck/carina/Tomato/Chrom_FeatureOverlap/ALLGENES.txt",
#quiet=TRUE, what="" )

K27posGenes=scan(
  "Z:/dep_coupland/grp_turck/carina/Tomato/Chrom_FeatureOverlap/H3K27me3+genes.txt",
  quiet=TRUE, what="")

Gene <- list(K27posGenes)
Genebackground<-c("K27+_Genes")
TargetGeneClusters=list.files(pattern="\\genes_GenesinClusters\\.txt$", recursive=TRUE)
#TargetCDSClusters=list.files(pattern="\\CDSinClusters\\.txt$", recursive=TRUE)

RY=c("cluster_1","cluster_2", "cluster_3", "cluster_4")
Telobox=c("cluster_1","cluster_2", "cluster_3", "cluster_4")
Telolike=c("cluster_1","cluster_2","cluster_3", "cluster_4")

#needed in loop: 
#unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t", col.names = c("CHR", "gene","cluster")) %>%
# rowwise() %>%
#mutate(gene=gsub("\\.[1-9]\\.[1-9]$","",gene))

#### Calculations ####

Background<-lengths(Gene)   
Organism <- c("Slycopersicum")

for (i in 1:length(TargetGeneClusters)) {
  
  unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=FALSE, col.names = c("Genome","gene", "cluster"))%>%
    rowwise() %>%
    mutate(gene=gsub("\\.[1-9]\\.[1-9]$","",gene)) 
  
  
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

####Medicago truncatula ####

setwd("Z:/dep_coupland/grp_turck/carina/Mtruncatula/Downloads")

geneID2GO <- readMappings("MtrunA17r5.0-ANR-EGN-r1.8.b2g.topgo")

setwd("Z:/dep_coupland/grp_turck/carina/Mtruncatula/GO-Analysis")

#Gene files for TopGO: scan Gene txt files. Create vector with p-values and assign gene names to vector members 
Genes=scan(
  "Mtruncatula_Genes.txt",
  quiet=TRUE,
  what=""
)

K27posGenes=scan(
  "Z:/dep_coupland/grp_turck/carina/Mtruncatula/Chrom_FeatureOverlap/Genes/H3K27me3K27+Genes.txt",
  quiet=TRUE,
  what=""
)

Gene <- list(K27posGenes)
Genebackground<-c("K27+_Genes")
TargetGeneClusters=list.files(pattern="\\GenesinClusters\\.txt$", recursive=TRUE)
#TargetCDSClusters=list.files(pattern="\\CDSinClusters\\.txt$", recursive=TRUE)

RY=c("cluster_1","cluster_2","cluster_3","cluster_4")
Telobox=c("cluster_1","cluster_2","cluster_3","cluster_4")
Telolike=c("cluster_1","cluster_2","cluster_3","cluster_4")

#### Calculations ####
Background<-lengths(Gene)   
Organism <- c("Mtruncatula")

for (i in 1:length(TargetGeneClusters)) {
  
  unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster")) 
  
  
  splitfilename=strsplit(TargetGeneClusters[[i]], '[-_|.]')[[1]]
  teloclustername=paste0(splitfilename[1],"_",splitfilename[2])
  
  if  ( splitfilename[1] == "RY" ) {
    TargetGenes_RY <- nrow(unfilteredTargetGenes %>% filter(cluster %in% RY)) 
    Background_RY <-nrow(unfilteredTargetGenes)
  } else if ( splitfilename[1]  ==  "Telobox" | splitfilename[1] == "telobox") {
    TargetGenes_TB <- nrow(unfilteredTargetGenes %>% filter(cluster %in% Telobox)) 
    Background_TB <- nrow(unfilteredTargetGenes)
  } else if ( splitfilename[1]  ==  "Telolike" | splitfilename[1] == "telolike" ) {
    TargetGenes_TL <- nrow(unfilteredTargetGenes %>% filter(cluster %in% Telolike)) 
    Background_TL <-nrow(unfilteredTargetGenes)
  }
}
newGennumbers<-data.frame(Organism,Background,TargetGenes_RY,TargetGenes_TB,TargetGenes_TL)

Genenumbers<- rbind(Genenumbers,newGennumbers)

###File####

setwd("~/Studium/Master/Master thesis/Projektmodul Turck/GO-Analyse")

write.xlsx(Genenumbers, file = "GO_Genenumbers.xlsx")