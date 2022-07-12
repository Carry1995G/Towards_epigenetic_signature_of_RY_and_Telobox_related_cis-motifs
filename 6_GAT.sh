#!/bin/bash

cd ~

#variables and paths
path1=$(pwd)
mkdir -p ${path1}/GAT
GAT=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/GAT
IDR=${path1}/IDR
TELOBED=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Motif_BED
EPIC=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Peaks
ANNOTATION=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Annotation_Files
DOWNLOADS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Downloads
BLACKLISTS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists

##Main Workspace##
bedtools intersect -v -a ${path1}/Extras/Mtruncatula_chromlen.bed -b ${BLACKLISTS}/Contigs_PlastidBL.bed >${GAT}/GAT_workspace.bed 

###Other workspaces, e.g. used in S. lycopersicum ###

#workspace minus Repeatregions + telobox rich regions (Blacklist) #
#gff2bed < ${DOWNLOADS}/*repeats*.gff > ${DOWNLOADS}/Aggressive_repeats.bed
#tail -n 12 ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_CHROMs.bed
#bedtools subtract -a GAT_workspace_CHROMs.bed -b ${DOWNLOADS}/Aggressive_repeats.bed > ${GAT}/GAT_workspace_norepeats.bed
#bedtools subtract -a ${GAT}/GAT_workspace_norepeats.bed -b ${BLACKLISTS}/Mtruncatula_telomerregions_BL.bed > ${GAT}/GAT_workspace_repeats+BL.bed
#bedtools subtract -a ${GAT}/GAT_workspace_CHROMs.bed -b ${ANNOTATION}/GFF3_intergenic.bed > ${GAT}/GAT_workspace_nointergenic.bed
#bedtools subtract -a ${GAT}/GAT_workspace_norepeats.bed -b ${ANNOTATION}/GFF3_intergenic.bed > ${GAT}/GAT_workspace_norpts+intergenic.bed
#bedtools subtract -a ${GAT}/GAT_workspace_repeats+BL.bed -b ${ANNOTATION}/GFF3_intergenic.bed > ${GAT}/GAT_workspace_norpts+intergenic+telorpts


###Prepare Annotations###
bedtools intersect -f 0.5 -a ${ANNOTATION}/GAT.bed -b ${GAT}/GAT_workspace.bed > ${GAT}/Annot.bed
awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $8, $8, $6}' ${GAT}/Annot.bed > ${GAT}/GAT.bed  #change columns into correct format for GAT
Annot=(${GAT}/GAT.bed)
bedtools intersect -f 0.5 -u -a  ${Annot} -b ${EPIC}/*intersect.bed | awk '{$4=$4 "_K27me3+"} {$5=$4} 1' OFS='\t'> ${GAT}/K27+_Annot_GAT.bed
bedtools intersect -f 0.5 -v -a ${Annot} -b  ${EPIC}/*intersect.bed | awk '{$4=$4 "_K27me3-"} {$5=$4} 1' OFS='\t'> ${GAT}/K27-_Annot_GAT.bed
cat ${GAT}/*Annot_GAT.bed | sort -k1,1 -k2,2n > ${GAT}/K27_all_GAT.bed


#Loop for Gat analysis to use each chromosome by itself#
#for z in $(seq 1 $(wc -l < ${GAT}/GAT_workspace.bed))
#do
#sed -n ${z}p ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_CHROM.bed

### GAT Analyse ####
rm -f ${GAT}/gat.log

cd ${TELOBED} 				#move to folder with Motif-BED files

for i in *.bed
do
gat-run.py --ignore-segment-tracks \
	--segments=$i \
	--annotations=${GAT}/K27_all_GAT.bed \
	--workspace=${GAT}/GAT_workspace.bed\
	--num-samples=1000 \
	--log=${GAT}/gat.log > ${GAT}/${i%_Mtruncatula.bed}gat.out
awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%_Mtruncatula.bed}gat.out >${GAT}/${i%.bed}table.txt
awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%_Mtruncatula.bed}gat.out >${GAT}/${i%.bed}header.txt
cat ${GAT}/${i%_Mtruncatula.bed}header.txt ${GAT}/${i%_Mtruncatula.bed}table.txt > ${GAT}/${i%_Mtruncatula.bed}GAT_result.txt
done

rm -f ${GAT}/*table.txt
rm -f  ${GAT}/*header.txt
rm -f ${GAT}/*gat.out
