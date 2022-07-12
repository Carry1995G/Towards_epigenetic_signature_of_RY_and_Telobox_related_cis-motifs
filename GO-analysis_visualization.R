##GO-Analysis Visualization for multiple organisms#

# Required packages
library(ggplot2)
library(xlsx)
library(dplyr)
library(pgirmess)
library(tidyr)
library(stringr)
library(scales)

# set the working directory where the tables to use are located
setwd("GO-Analyse")
motif=c("Telobox", "Telolike", "RY")
background=c("K27[+][_]Genes")
organism="xxx"            #change accordingly

## Plotting ##
for (m in motif) {
  GOfiles=list.files(pattern=sprintf("%s_%s_%s.xlsx",m,organism,background))
  
  GO_gp1<-read.xlsx(GOfiles[1],sheetIndex=1 ,header=T,stringsAsFactors = T)  %>% slice_head(n = 5)
  GO_gp2<-read.xlsx(GOfiles[2],sheetIndex=1 ,header=T,stringsAsFactors = T)  %>% slice_head(n = 5)
  GO_gp3<-read.xlsx(GOfiles[3],sheetIndex=1 ,header=T,stringsAsFactors = T)  %>% slice_head(n = 5)
  
  splitname=strsplit(GOfiles[[1]], '[-_]')[[1]]
  GOprocess=splitname[1]
  GO_gp1 <- GO_gp1 %>% mutate(Process =strsplit(GOfiles[[1]], '[-_]')[[1]][1] )
  GO_gp2 <- GO_gp2 %>% mutate(Process =strsplit(GOfiles[[2]], '[-_]')[[1]][1] )
  GO_gp3 <- GO_gp3 %>% mutate(Process =strsplit(GOfiles[[3]], '[-_]')[[1]][1] )
  
  GO_gp=rbind(GO_gp1,GO_gp2,GO_gp3) %>% mutate(Fold_change=(Significant/Expected)-1) %>% 
    filter(pvalue<0.01) %>%
    mutate(Annotated=as.numeric(Annotated),# Transform the column 'Gene_number' into a numeric variable
           GO.ID = paste("(",GO.ID,")",sep="")) %>%
    mutate(Term=as.character(Term),
           Term=str_trunc(Term,23,"right")) %>%#include () around GO-Terms
    unite(GO_Term, Term, GO.ID, sep= " ") %>% #combine GO-Terms and Processes
    mutate (GO_Term=as.factor(GO_Term), # Transform the column 'GO_biological_process' into factors
            '|-log10(pvalue)|'=-(log10(pvalue)))  # Transform padj values by -log10('padj')

  # Draw the plot in facets by group with ggplot2
  # to represent -log10(FDR), Number of genes and 
  # Fold enrichment of each GO biological process per group (Figure 3)
  #-------------------------------------------------------------------
  ggplot(GO_gp, aes(y = GO_Term, x= Fold_change)) +
    geom_point(data=GO_gp,aes(y=GO_Term, x= Fold_change,size = Significant, colour =`|-log10(pvalue)|` ), alpha=.7)+
    scale_color_gradient(low="green",high="red",limits=c(0, NA),breaks=breaks_extended(n=5))+ 
    scale_size_continuous(breaks=breaks_extended(n=4))+
    theme_bw()+
    theme(axis.ticks.length=unit(-0.1, "cm"),
          axis.text.x = element_text(margin=margin(5,5,0,5,"pt")),
          axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
          axis.text = element_text(color = "black"),
          panel.grid.minor = element_blank(),
          legend.title.align=0.5, legend.title = element_text(size=12),
          text=element_text(size = 12))+
    ylab("GO Term")+
    xlab("Fold change")+
    ggtitle(paste(m,"motif"))+theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE),breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+#breaks=breaks_extended(n = 3)+
    labs(color="-log10\n(pvalue)", size="Number\nof genes") +
    facet_wrap(~Process, strip.position = "bottom",ncol=4)
  
  ggsave (filename=paste(m,"_K27genes_pvalueGOAnalysis.png"), width=5, height=4)
}

remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)
