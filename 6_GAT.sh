#!/bin/bash

cd ~

#variables and paths
path1=$(pwd)
mkdir -p ${path1}/GAT
GAT=/GAT
TELOBED=/Motif_BED
EPIC=/Peaks
ANNOTATION=/Annotation_Files
DOWNLOADS=/Downloads
BLACKLISTS=/Blacklists

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

### GAT Analysis ####
rm -f ${GAT}/gat.log

cd ${TELOBED} 				#move to folder with Motif-BED files

for i in *.bed
do
gat-run.py --ignore-segment-tracks \
	--segments=$i \
	--annotations=${GAT}/K27_all_GAT.bed \
	--workspace=${GAT}/GAT_workspace.bed\
	--num-samples=1000 \
	--log=${GAT}/gat.log > ${GAT}/${i%.bed}gat.out
awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%.bed}gat.out >${GAT}/${i%.bed}table.txt
awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%.bed}gat.out >${GAT}/${i%.bed}header.txt
cat ${GAT}/${i%.bed}header.txt ${GAT}/${i%.bed}table.txt > ${GAT}/${i%.bed}GAT_result.txt
done

rm -f ${GAT}/*table.txt ${GAT}/*header.txt ${GAT}/*gat.out

#GAT Analysis used for tomato#
#Loop to use each chromosome by itself. Used for S. lypcopersicum#
for z in $(seq 1 $(wc -l < ${GAT}/GAT_workspace.bed))   #for z in length of workspace file
do
sed -n ${z}p ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_CHROM.bed #choose line nr z

for i in *.bed
do
gat-run.py --ignore-segment-tracks \
	--segments=$i \
	--annotations=${GAT}/K27_all_GAT.bed \
	--workspace=${GAT}/${GAT}/GAT_workspace_CHROM.bed\
	--num-samples=1000 \
	--log=${GAT}/gat.log > ${GAT}/${i%.bed}gat.out
awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%.bed}gat.out >${GAT}/${i%.bed}table.txt
awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%.bed}gat.out >${GAT}/${i%.bed}header.txt
cat ${GAT}/${i%.bed}header.txt ${GAT}/${i%.bed}table.txt > ${GAT}/${i%.bed}GAT_result.txt
done
done 

rm -f ${GAT}/*table.txt ${GAT}/*header.txt ${GAT}/*gat.out
