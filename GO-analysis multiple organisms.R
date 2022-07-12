####GO-Analysis for multiple organisms#######

#Preparation
library(topGO)
library(pgirmess)
library(xlsx)
library(dplyr)
library(tidyr)
setwd("~/Studium/Master/Master thesis/Projektmodul Turck/GO-Analyse")

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#############################TopGO#########################

###A thaliana####
setwd("Z:/dep_coupland/grp_turck/carina/Athaliana/GO-Analysis")

geneID2GO <- (readMappings("GO_Map"))

#Genes=scan(
 # "genes.txt",
  #quiet=TRUE,
  #what=""
#)

K27posGenes<-scan(
  "Z:/dep_coupland/grp_turck/carina/Athaliana/Chrom_FeatureOverlap/H3K27me3+_genes.txt",
  quiet=TRUE,
  what=""
)

RY=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telobox=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telolike=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")

Gene <- list(K27posGenes)
Genebackground<-c("K27+_Genes")
TargetGeneClusters=list.files(pattern="\\GenesinClusters\\.txt$", recursive=TRUE)

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


###B distachyon####

setwd("Z:/dep_coupland/grp_turck/carina/Bdistachyon/GO-Analysis")

geneID2GO <- (readMappings("GO_Map"))

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
####Solanum lycopersicum####

setwd("Z:/dep_coupland/grp_turck/carina/Tomato/GO-Analysis")

geneID2GO <- readMappings("Z:/dep_coupland/grp_turck/carina/Tomato/Downloads/ITAG4.0_goterms.txt")


#Gene files for TopGO: scan Gene txt files. Create vector with p-values and assign gene names to vector members 
#Genes=scan(
 # "Z:/dep_coupland/grp_turck/carina/Tomato/Chrom_FeatureOverlap/ALLGENES.txt",
  #quiet=TRUE, what="" )

K27posGenes=scan(
  "Z:/dep_coupland/grp_turck/carina/Tomato/Chrom_FeatureOverlap/H3K27me3+genes.txt",
  quiet=TRUE, what="")

#Gene <- list(K27posGenes)
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

for (i in 1:length(TargetGeneClusters)) {
    
 unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster")) 
  
  #Marchantia:
  #unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster"))  %>% rowwise() %>% #, col.names = c("gene", "cluster")) 
   #mutate(gene=gsub(".*\\M","M",gene), 
    #     gene=ifelse(gene %in% translationshort$V2[duplicated(translationshort$V2)],"DOUBLE",
     #               ifelse(gene %in% unique(translationshort$V2), translationshort$V1[which(translationshort$V2==gene)],
      #                    gene)))
 
 #Solanum lycopersicum
 #unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t", col.names = c("CHR", "gene","cluster")) %>%
  # rowwise() %>%
   #mutate(gene=gsub("\\.[1-9]\\.[1-9]$","",gene))
  
  splitfilename=strsplit(TargetGeneClusters[[i]], '[-_|.]')[[1]]
  teloclustername=paste0(splitfilename[1],"_",splitfilename[2])
  
  if  ( splitfilename[1] == "RY" ) {
   TargetGenes <-  unfilteredTargetGenes %>% filter(cluster %in% RY) 
  } else if ( splitfilename[1]  ==  "Telobox" | splitfilename[1] == "telobox") {
    TargetGenes <- unfilteredTargetGenes %>% filter(cluster %in% Telobox) 
  } else if ( splitfilename[1]  ==  "Telolike" | splitfilename[1] == "telolike" ) {
    TargetGenes <- unfilteredTargetGenes %>% filter(cluster %in% Telolike) 
  }
  
  ###Gene files with p-values
  #pvalue_targetGenes = rep(c(0.001),times=length(TargetGenes))
  #names(pvalue_targetGenes)=c(TargetGenes)
  #pvalue_allGenes = rep(c(1),times=19213)
  #names(pvalue_allGenes)=c(Genes)   ###assign names to vector members
  
  ####predefined list of interesting genes
  for (g in 1) {
  geneNames <- Gene[[g]] #names(geneID2GO)
  myInterestingGenes <- c(TargetGenes$gene)
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  
  ##store all data in sampleGOdata to faciliate access to identifiers, annotation 
  ##and basic data statistics
  
  sampleMFGOdata <- new("topGOdata", ontology = "MF",
                        allGenes = geneList, annot = annFUN.gene2GO, 
                        gene2GO=geneID2GO)
  
  sampleCCGOdata <- new("topGOdata", ontology = "CC",
                        allGenes = geneList, annot = annFUN.gene2GO, 
                        gene2GO=geneID2GO)
  
  sampleBPGOdata <- new("topGOdata", ontology = "BP",
                        allGenes = geneList, annot = annFUN.gene2GO, 
                        gene2GO=geneID2GO)
  
  ###########Fisher-Test
  #Selecting 'algorithm=classic' means that the GO hierarchy isn't taken into account, 
  #so each GO term is tested independently (over-representation enrichment). 
  #limitation: all genes annotated to a GO term will be automatically annotated to its parents as well
  #therefore a GO term might look enriched just because its children are enriched. 
  #Thus, it is important that GO hierarchy is taken into account (conditional enrichment) to avoid redundancy.
  resultFisMF <- runTest(sampleMFGOdata, algorithm='weight01', statistic = "fisher")
  resultFisCC <- runTest(sampleCCGOdata, algorithm='weight01', statistic = "fisher")
  resultFisBP <- runTest(sampleBPGOdata, algorithm='weight01', statistic = "fisher")
  
  ######## CreateTable######
  
  concepts <- list("BP"=resultFisBP, "CC"=resultFisCC, "MF"=resultFisMF)
  
  for (a in 1:length(concepts)) { 
    
    Bioconcept = concepts[[a]]
    ABR = names(concepts[a])
    
    if  ( ABR ==  "MF" ) {
      GOdata = sampleMFGOdata
    } else if ( ABR ==  "CC" ) {
      GOdata = sampleCCGOdata
    } else if ( ABR ==  "BP" ) {
      GOdata = sampleBPGOdata
    }
    
    ###Calculate FDR 
    allGO = usedGO(object = GOdata)
    score <- data.frame(GO.ID=names(score(Bioconcept)), 
                        pvalue=score(Bioconcept), row.names=NULL)

    allRes <-GenTable(GOdata, weightFisher=Bioconcept, orderBy='weightFisher', topNodes = length(allGO))
    allRes <- left_join(allRes, score, by = "GO.ID")  
    p.adj=p.adjust(allRes$pvalue, method="BH",n=length(allRes$pvalue))
    allRes$p.adj=p.adj

    #get list of significant GO after multiple testing correction 0.05)
    results.BH=allRes[which(allRes$p.adj <= 0.05),]
    
    #get list of significant GO before multiple testing correction
    results.p <- allRes[which(allRes$pvalue<=0.05),]

    #create columns for table
    GO_biological_process=paste0(results.BH$Term,"(",results.BH$GO.ID,")")
    GO_biological_process=gsub(" ", "_", GO_biological_process, fixed = TRUE)
    
    Fold_enrichment = results.BH$Significant/results.BH$Expected
    FDR=results.BH$p.adj
    Gene_number=results.BH$Significant
    Group=rep(teloclustername,times=nrow(results.BH))
    ##export table
    write.xlsx(allRes, file = paste0(ABR,"_",teloclustername,"_",Genebackground[g],".xlsx"), sheetName="Full_Results")
    print ("write xlsx done")
    sink (file = paste0(ABR,"_",teloclustername,"_",Genebackground[g],".txt"))
    print(Bioconcept)
    sink()
    
    if ( nrow(results.p) == 0) { 
      if (nrow(results.BH) == 0) { 
        next 
      } else { 
        
        table=data.frame(GO_biological_process,	Gene_number, Fold_enrichment,	FDR)
        write.xlsx(results.BH, file = paste0(ABR,"_",teloclustername,"_",Genebackground[g],".xlsx"), sheetName="BH_results", append=TRUE)
        write.xlsx(table, file = paste0(ABR,"_",teloclustername,"_",Genebackground[g],".xlsx"),  sheetName="FDR<0.01", append=TRUE )
      }
    } else { 
      write.xlsx(results.p, file = paste0(ABR,"_",teloclustername,"_",Genebackground[g],".xlsx"), sheetName="pvalue_results",append=TRUE)
      if (nrow(results.BH) == 0) { 
        next 
      } else { 
        
        table=data.frame(GO_biological_process,	Gene_number, Fold_enrichment,	FDR)
        write.xlsx(results.BH, file = paste0(ABR,"_",teloclustername,"_",Genebackground[g],".xlsx"), sheetName="BH_results", append=TRUE)
        write.xlsx(table, file = paste0(ABR,"_",teloclustername,"_",Genebackground[g],".xlsx"),  sheetName="FDR<0.01", append=TRUE )
        
      }
    
      ####export group table
      grouptable=data.frame(Group,GO_biological_process,	Gene_number,	Fold_enrichment,	FDR)
      
      # assign flexible variable with grouptable within loop
      assign(paste0(ABR,"_group_",teloclustername), grouptable)
    }
  }
  }
}











####Grouping - not necessary ####
###just if analysis makes difference between clusters ###

group_MF_telolike = rbind(MF_group_telolike_cluster1, MF_group_telolike_cluster2, MF_group_telolike_cluster4)
group_CC_telolike =rbind(CC_group_telolike_cluster1,CC_group_telolike_cluster2,CC_group_telolike_cluster4)
group_BP_telolike=rbind(BP_group_telolike_cluster1, BP_group_telolike_cluster4)
group_MF_telobox = rbind(MF_group_telobox_cluster1, MF_group_telobox_cluster2, MF_group_telobox_cluster4)
group_CC_telobox =rbind (CC_group_telobox_cluster1,CC_group_telobox_cluster2,CC_group_telobox_cluster4)
group_BP_telobox = rbind(BP_group_telobox_cluster1, BP_group_telobox_cluster4)
group_MF_RY<-rbind(BP_group_telobox_cluster1, BP_group_telobox_cluster4)
group_CC_RY<-rbind(BP_group_telobox_cluster1, BP_group_telobox_cluster4)
group_BP_RY<-rbind(BP_group_telobox_cluster1, BP_group_telobox_cluster4)


groupfiles = list("group_telobox_MF"=group_MF_telobox, 
                  "group_telobox_CC"=group_CC_telobox,
                  "group_telobox_BP"=group_BP_telobox,
                  "group_telolike_MF"=group_MF_telolike,
                  "group_telolike_CC"=group_CC_telolike,
                  "group_telolike_BP"=group_BP_telolike, 
                  "group_RY_MF"=group_MF_RY,
                  "group_RY_CC"=group_CC_RY,
                  "group_RY_BP"=group_BP_RY)

for (a in 1:length(groupfiles)) {
  fullgrouptable = groupfiles[[a]]
  filename = names(groupfiles[a])
  
  write.delim(fullgrouptable, file = paste0(filename,".txt"), row.names = FALSE, quote = FALSE, sep = "\t" )
}