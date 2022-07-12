##GO-Analysis Visualization for multiple organisms#

# Requires the package 'ggplot2' (needs to be installed first)
# Load the ggplot2 package
library(ggplot2)
library(xlsx)
library(dplyr)
library(pgirmess)
library(tidyr)
library(stringr)
library(scales)


# set the working directory where the tables to use are located
setwd("~/Studium/Master/Projektmodul Turck/GO-Analyse")


###Athaliana####
setwd("Z:/dep_coupland/grp_turck/carina/Athaliana/GO-Analysis")

motif=c("Telobox", "Telolike", "RY")
background=c("K27[+][_]Genes")
organism="Athaliana"

#### Plotting ####
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
  # Import the table containing the enriched GO terms by groups
  #GO_gp <- read.table(GOfiles[1],header=T,stringsAsFactors = T, fill = TRUE )
  
  # List objects and their structure contained in the dataframe 'GO_gp'
  ##ls.str(GO_gp)
  
  # Change factor order
  #GO_gp$Group<- factor(GO_gp$Group,levels = c("WT_up","WT_down","mutant_up","mutant_down"))
  # GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))
  
  # Create a vector with new names for groups to use in the plot
  # Replace the terms by your own (\n allow to start a new line)
  #group.labs <- c(`WT_up` = "WT up-\nregulated",
  #               `WT_down` = "WT down-\nregulated",
  #              `mutant_up` = "Mutant up-\nregulated",
  #             `mutant_down` = "Mutant down-\nregulated")
  
  # Draw the plot in facets by group with ggplot2
  # to represent -log10(FDR), Number of genes and 
  # Fold enrichment of each GO biological process per group (Figure 3)
  #-------------------------------------------------------------------
  ggplot(GO_gp, aes(y = GO_Term, x= Fold_change)) +
    # geom_vline(xintercept = 1, linetype="dashed", 
    #           color = "azure4", size=.5)+
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
    facet_wrap(~Process, strip.position = "bottom",ncol=4)#labeller=as_labeller(group.labs))+#after "ncol=", specify the number of groups you have
  #guides(y = guide_legend(order=2),
  #  colour = guide_colourbar(order=1))
  
  ggsave (filename=paste(m,"_K27genes_pvalueGOAnalysis.png"), width=5, height=4)
}

remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

####Marchantia polymorpha####
remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

setwd("Z:/dep_coupland/grp_turck/carina/Moss-Liverworts/Marchantia-chromatin/GO-Analyse")

motif=c("RY")
background=c("K[0-9]{2}[+][_]Genes")
organism="Mpolymorpha"


#### Plotting ####
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
  
  GO_gp=rbind(GO_gp1,GO_gp2,GO_gp3) %>% mutate(Fold_change=Significant/Expected-1) %>%
    filter(pvalue<0.01) %>%
    mutate(Annotated=as.numeric(Annotated),# Transform the column 'Gene_number' into a numeric variable
           GO.ID = paste("(",GO.ID,")",sep="")) %>%
    mutate(Term=as.character(Term),
           Term=str_trunc(Term,23,"right")) %>%#include () around GO-Terms
    unite(GO_Term, Term, GO.ID, sep= " ") %>% #combine GO-Terms and Processes
    mutate (GO_Term=as.factor(GO_Term), # Transform the column 'GO_biological_process' into factors
            '|-log10(pvalue)|'=-(log10(pvalue)))  # Transform padj values by -log10('padj')
  # Import the table containing the enriched GO terms by groups
  #GO_gp <- read.table(GOfiles[1],header=T,stringsAsFactors = T, fill = TRUE )
  
  # List objects and their structure contained in the dataframe 'GO_gp'
  #ls.str(GO_gp)
  
  # Change factor order
  #GO_gp$Group<- factor(GO_gp$Group,levels = c("WT_up","WT_down","mutant_up","mutant_down"))
  # GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))
  
  # Create a vector with new names for groups to use in the plot
  # Replace the terms by your own (\n allow to start a new line)
  #group.labs <- c(`WT_up` = "WT up-\nregulated",
  #               `WT_down` = "WT down-\nregulated",
  #              `mutant_up` = "Mutant up-\nregulated",
  #             `mutant_down` = "Mutant down-\nregulated")
  
  # Draw the plot in facets by group with ggplot2
  # to represent -log10(FDR), Number of genes and 
  # Fold enrichment of each GO biological process per group (Figure 3)
  #-------------------------------------------------------------------
  ggplot(GO_gp, aes(y = GO_Term, x= Fold_change)) +
    # geom_vline(xintercept = 1, linetype="dashed", 
    #           color = "azure4", size=.5)+
    geom_point(data=GO_gp,aes(y=GO_Term, x= Fold_change,size = Significant, colour =`|-log10(pvalue)|` ), alpha=.7)+
    
    scale_color_gradient(low="green",high="red",limits=c(0, NA),breaks=breaks_extended(n=5))+ scale_size_continuous(breaks=breaks_extended(n=4))+     #coord_flip()+
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
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
    labs(color="-log10\n(pvalue)", size="Number\nof genes") +
    facet_wrap(~Process, strip.position = "bottom",ncol=4)#labeller=as_labeller(group.labs))+#after "ncol=", specify the number of groups you have
  #guides(y = guide_legend(order=2),
  #  colour = guide_colourbar(order=1))
  
  ggsave (filename=paste(m,"_K27genes_pvalueGOAnalysis.png"), width=5, height=4)
}

remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

###Bdistachyon####
remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

setwd("Z:/dep_coupland/grp_turck/carina/Bdistachyon/GO-Analysis")

motif=c("Telobox", "RY")
background=c("K27[+][_]Genes")
organism="Bdistachyon"

#### Plotting ####
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
  
  GO_gp=rbind(GO_gp1,GO_gp2,GO_gp3) %>% mutate(Fold_change=Significant/Expected-1) %>% 
    filter(pvalue<0.01) %>%
    mutate(Annotated=as.numeric(Annotated),# Transform the column 'Gene_number' into a numeric variable
           GO.ID = paste("(",GO.ID,")",sep="")) %>%
    mutate(Term=as.character(Term),
           Term=str_trunc(Term,23,"right")) %>%#include () around GO-Terms
    unite(GO_Term, Term, GO.ID, sep= " ") %>% #combine GO-Terms and Processes
    mutate (GO_Term=as.factor(GO_Term), # Transform the column 'GO_biological_process' into factors
            '|-log10(pvalue)|'=-(log10(pvalue)))  # Transform padj values by -log10('padj')
  # Import the table containing the enriched GO terms by groups
  #GO_gp <- read.table(GOfiles[1],header=T,stringsAsFactors = T, fill = TRUE )
  
  # List objects and their structure contained in the dataframe 'GO_gp'
  #ls.str(GO_gp)
  
  # Change factor order
  #GO_gp$Group<- factor(GO_gp$Group,levels = c("WT_up","WT_down","mutant_up","mutant_down"))
  # GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))
  
  # Create a vector with new names for groups to use in the plot
  # Replace the terms by your own (\n allow to start a new line)
  #group.labs <- c(`WT_up` = "WT up-\nregulated",
  #               `WT_down` = "WT down-\nregulated",
  #              `mutant_up` = "Mutant up-\nregulated",
  #             `mutant_down` = "Mutant down-\nregulated")
  
  # Draw the plot in facets by group with ggplot2
  # to represent -log10(FDR), Number of genes and 
  # Fold enrichment of each GO biological process per group (Figure 3)
  #-------------------------------------------------------------------
  ggplot(GO_gp, aes(y = GO_Term, x= Fold_change)) +
    # geom_vline(xintercept = 1, linetype="dashed", 
    #           color = "azure4", size=.5)+
    geom_point(data=GO_gp,aes(y=GO_Term, x= Fold_change,size = Significant, colour =`|-log10(pvalue)|` ), alpha=.7)+
    
    scale_color_gradient(low="green",high="red",limits=c(0, NA),breaks=breaks_extended(n=5))+ scale_size_continuous(breaks=breaks_extended(n=4))+     #coord_flip()+
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
    scale_x_continuous(breaks = breaks_extended(n=3), guide = guide_axis(check.overlap = TRUE)) +# function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
    labs(color="-log10\n(pvalue)", size="Number\nof genes") +
    facet_wrap(~Process, strip.position = "bottom",ncol=4)#labeller=as_labeller(group.labs))+#after "ncol=", specify the number of groups you have
  #guides(y = guide_legend(order=2),
  #  colour = guide_colourbar(order=1))
  
  ggsave (filename=paste(m,"_K27genes_pvalueGOAnalysis.png"), width=5, height=4)
}

remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

###Tomato######
remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

setwd("Z:/dep_coupland/grp_turck/carina/Tomato/GO-Analysis")
motif=c("Telobox", "Telolike", "RY")
background=c("K[0-9]{2}[+][_]Genes")
organism="Tomato"

#### Plotting ####
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
  
  GO_gp=rbind(GO_gp1,GO_gp2,GO_gp3) %>% mutate(Fold_change=Significant/Expected-1) %>% 
    filter(pvalue<0.01) %>%
    mutate(Annotated=as.numeric(Annotated),# Transform the column 'Gene_number' into a numeric variable
           GO.ID = paste("(",GO.ID,")",sep="")) %>%
    mutate(Term=as.character(Term),
           Term=str_trunc(Term,23,"right")) %>%#include () around GO-Terms
    unite(GO_Term, Term, GO.ID, sep= " ") %>% #combine GO-Terms and Processes
    mutate (GO_Term=as.factor(GO_Term), # Transform the column 'GO_biological_process' into factors
            '|-log10(pvalue)|'=-(log10(pvalue)))  # Transform padj values by -log10('padj')
  # Import the table containing the enriched GO terms by groups
  #GO_gp <- read.table(GOfiles[1],header=T,stringsAsFactors = T, fill = TRUE )
  
  # List objects and their structure contained in the dataframe 'GO_gp'
  #ls.str(GO_gp)
  
  # Change factor order
  #GO_gp$Group<- factor(GO_gp$Group,levels = c("WT_up","WT_down","mutant_up","mutant_down"))
  # GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))
  
  # Create a vector with new names for groups to use in the plot
  # Replace the terms by your own (\n allow to start a new line)
  #group.labs <- c(`WT_up` = "WT up-\nregulated",
  #               `WT_down` = "WT down-\nregulated",
  #              `mutant_up` = "Mutant up-\nregulated",
  #             `mutant_down` = "Mutant down-\nregulated")
  
  # Draw the plot in facets by group with ggplot2
  # to represent -log10(FDR), Number of genes and 
  # Fold enrichment of each GO biological process per group (Figure 3)
  #-------------------------------------------------------------------
  ggplot(GO_gp, aes(y = GO_Term, x= Fold_change)) +
    # geom_vline(xintercept = 1, linetype="dashed", 
    #           color = "azure4", size=.5)+
    geom_point(data=GO_gp,aes(y=GO_Term, x= Fold_change,size = Significant, colour =`|-log10(pvalue)|` ), alpha=.7)+
    
    scale_color_gradient(low="green",high="red",limits=c(0, NA),breaks=breaks_extended(n=5))+ scale_size_continuous(breaks=breaks_extended(n=4))+     #coord_flip()+
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
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
    labs(color="-log10\n(pvalue)", size="Number\nof genes") +
    facet_wrap(~Process, strip.position = "bottom",ncol=4)#labeller=as_labeller(group.labs))+#after "ncol=", specify the number of groups you have
  #guides(y = guide_legend(order=2),
  #  colour = guide_colourbar(order=1))
  
  ggsave (filename=paste(m,"_K27genes_pvalueGOAnalysis.png"), width=5, height=4)
}

remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

####Medicago truncatula####
remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)

setwd("Z:/dep_coupland/grp_turck/carina/Mtruncatula/GO-Analysis")

motif=c("Telobox", "Telolike", "RY")
background=c("K[0-9]{2}[+][_]Genes")
organism="Mtruncatula"


#### Plotting ####
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
  
  GO_gp=rbind(GO_gp1,GO_gp2,GO_gp3) %>% mutate(Fold_change=Significant/Expected-1) %>% 
    filter(pvalue<0.01) %>%
    mutate(Annotated=as.numeric(Annotated),# Transform the column 'Gene_number' into a numeric variable
           GO.ID = paste("(",GO.ID,")",sep="")) %>%
    mutate(Term=as.character(Term),
           Term=str_trunc(Term,23,"right")) %>%#include () around GO-Terms
    unite(GO_Term, Term, GO.ID, sep= " ") %>% #combine GO-Terms and Processes
    mutate (GO_Term=as.factor(GO_Term), # Transform the column 'GO_biological_process' into factors
            '|-log10(pvalue)|'=-(log10(pvalue)))  # Transform padj values by -log10('padj')
  # Import the table containing the enriched GO terms by groups
  #GO_gp <- read.table(GOfiles[1],header=T,stringsAsFactors = T, fill = TRUE )
  
  # List objects and their structure contained in the dataframe 'GO_gp'
  #ls.str(GO_gp)
  
  # Change factor order
  #GO_gp$Group<- factor(GO_gp$Group,levels = c("WT_up","WT_down","mutant_up","mutant_down"))
  # GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=rev(levels(GO_gp$GO_biological_process)))
  
  # Create a vector with new names for groups to use in the plot
  # Replace the terms by your own (\n allow to start a new line)
  #group.labs <- c(`WT_up` = "WT up-\nregulated",
  #               `WT_down` = "WT down-\nregulated",
  #              `mutant_up` = "Mutant up-\nregulated",
  #             `mutant_down` = "Mutant down-\nregulated")
  
  # Draw the plot in facets by group with ggplot2
  # to represent -log10(FDR), Number of genes and 
  # Fold enrichment of each GO biological process per group (Figure 3)
  #-------------------------------------------------------------------
  ggplot(GO_gp, aes(y = GO_Term, x= Fold_change)) +
    # geom_vline(xintercept = 1, linetype="dashed", 
    #           color = "azure4", size=.5)+
    geom_point(data=GO_gp,aes(y=GO_Term, x= Fold_change,size = Significant, colour =`|-log10(pvalue)|` ), alpha=.7)+
    
    scale_color_gradient(low="green",high="red",limits=c(0, NA),breaks=breaks_extended(n=5))+ scale_size_continuous(breaks=breaks_extended(n=4))+     #coord_flip()+
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
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
    labs(color="-log10\n(pvalue)", size="Number\nof genes") +
    facet_wrap(~Process, strip.position = "bottom",ncol=4)#labeller=as_labeller(group.labs))+#after "ncol=", specify the number of groups you have
  #guides(y = guide_legend(order=2),
  #  colour = guide_colourbar(order=1))
  
  ggsave (filename=paste(m,"_K27genes_pvalueGOAnalysis.png"), width=5, height=4)
}

remove(GO_gp,GO_gp1,GO_gp2,GO_gp3)
