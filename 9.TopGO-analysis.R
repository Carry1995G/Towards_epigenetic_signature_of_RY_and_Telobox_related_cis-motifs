####TopGO-Analysis ####

#Preparation
library(topGO)
library(pgirmess)
library(xlsx)
library(dplyr)
library(tidyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

setwd("~/GO-Analysis")
geneID2GO <- (readMappings("GO_Map"))

K27posGenes<-scan(
  "Chrom_FeatureOverlap/H3K27me3+_genes.txt",
  quiet=TRUE,
  what=""
)

RY=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")              #Only contains clusters with high coverage
Telobox=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")
Telolike=c("cluster_1", "cluster_2", "cluster_3", "cluster_4")

Gene <- list(K27posGenes)
Genebackground<-c("K27+_Genes")
TargetGeneClusters=list.files(pattern="\\GenesinClusters\\.txt$", recursive=TRUE)   ##filename used from 8. GO_Preparation

#### Calculations ####

for (i in 1:length(TargetGeneClusters)) {
    
 unfilteredTargetGenes <- read.table(TargetGeneClusters[[i]],sep="\t",header=TRUE, col.names = c("gene", "cluster"))   ##file used from 8. GO_Preparation
  splitfilename=strsplit(TargetGeneClusters[[i]], '[-_|.]')[[1]]
  teloclustername=paste0(splitfilename[1],"_",splitfilename[2])
  
  if  ( splitfilename[1] == "RY" ) {
   TargetGenes <-  unfilteredTargetGenes %>% filter(cluster %in% RY) 
  } else if ( splitfilename[1]  ==  "Telobox" | splitfilename[1] == "telobox") {
    TargetGenes <- unfilteredTargetGenes %>% filter(cluster %in% Telobox) 
  } else if ( splitfilename[1]  ==  "Telolike" | splitfilename[1] == "telolike" ) {
    TargetGenes <- unfilteredTargetGenes %>% filter(cluster %in% Telolike) 
  }
  
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
    }
  }
  }
}
