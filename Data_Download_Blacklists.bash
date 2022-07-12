#!/bin/bash

###Author: Carina Gude
##Download genomic and chromatin data - 


#new folders and path links
mkdir -p Downloads
DOWNLOADS=/Downloads
cd ${DOWNLOADS}

####Obtain genomic data####

wget LINK

#if files are stored as .gz, unzip with:
gzip *.gz


#if files are stored as .zip files, unzip with 
unzip *.zip


###Obtain raw SRA data### 


###Generate Blacklist##
#depending on the genome, Contigs, Mitochondria and Plastids are stored in a blacklist to remove them from analysis. #

mkdir -p Blacklists
mkdir -p Extras
EXTRAS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Extras
BLACKLISTS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists

awk -F"\t" 'BEGIN {OFS="\t"} {print $0}' ${EXTRAS}/Chromlen.bed > ${Blacklists}/Contigs_PlastidBL.bed #BED files containing all chromosomes, contigs etc. 
sed -i "/Chromosome[1-x]/d" ${Blacklists}/Contigs_PlastidBL.bed   #remove chromosomes 1 to x. X=chromosome number. If chromosomes are named differently, a different pattern needs to be chosen
