#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/5_Blacklists.sh
cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

mkdir -p Blacklists
mkdir -p Plots
Extras=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Extras
Blacklists=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists
Telobed=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Motif_BED

####Write blacklist
cd ${Blacklists}

#Remove Chloroplast and Mitochondria from analysis: NOT NECESSARY FOR ITAG4.0
awk -F"\t" 'BEGIN {OFS="\t"} {print $0}' ${Extras}/Mtruncatula_chromlen.bed > ${Blacklists}/Contigs_PlastidBL.bed
sed -i "/MtrunA17Chr[1-8]/d" ${Blacklists}/Contigs_PlastidBL.bed

##Remove telomer regions
#bedtools makewindows -g ${Extras}/Mtruncatula_chromlen.txt -w 200 > Mtruncatula200bpwindows.bed
#bedtools coverage -a Mtruncatula200bpwindows.bed -b ${Telobed}/Telobox_Mtruncatula.bed -counts |  sort -k1,1 -k2,2n > Teloboxhitsper200bp.bed
#awk '$4>3' Teloboxhitsper200bp.bed > Mtruncatula_telomerregions_BL.bed

#Merge blacklists: not necessary for ITAG4.0
#cat Mtruncatula_telomerregions_BL.bed Min_Blacklist.bed >> merged_Blacklists.bed
#sort -k1,1 -k2,2n merged_Blacklists.bed > sorted_Blacklists.bed
#cut -f1-3 sorted_Blacklists.bed >> Telomer+Min_Blacklist.bed


#bsub -q bigmem -R "rusage[mem=20000]" -M 40000 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/11_GAT*.sh