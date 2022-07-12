#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/10_Replicate_Overlap.sh


###Description###
#Used instead of IDR #
#Careful: Bad replicate could result in less peaks for analysis #

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/
mkdir Peaks

#SRR14525744 - WT
#SRR14525745 - WT 
#SRR14525748 - WT H3K27me3
#SRR14525749 - WT H3K27me3

#cd BWA-align/BED-files
#mv SRR14525744.bed Control-R1
#mv SRR14525745.bed Control-R2
#mv SRR14525748.bed H3K27me3-R1
#mv SRR14525749.bed H3K27me3-R2

 
#variables and paths
path1=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula
path2=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Downloads/Pecrix*
EXTRAS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Extras
Blacklists=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists
PEAKS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Peaks

H3K27ME3=${path2}/H3K27me3_RNC1-R

mkdir -p ${path1}/IDR

#used:
awk -F"\t" -v OFS="\t" '{print $1,0,$2}' ${EXTRAS}/Mtruncatula_chromlen.txt >${EXTRAS}/Mtruncatula_chromlen.bed

cd ${path2}

bedtools intersect -a ${H3K27ME3}1/*islands-summary-FDR* -b ${H3K27ME3}2/*islands-summary-FDR* | bedtools intersect -wa -a - -b ${EXTRAS}/Mtruncatula_chromlen.bed |sort -k1,1 -k2,2n > ${PEAKS}/H3K27me3_SICER_intersect.bed