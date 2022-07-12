#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/12_GAT.sh
#bsub -q bigmem -R "rusage[mem=20000]" -M 40000 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/12_GAT.sh

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

####################################################################
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

echo "prepare workspace"
##Main Workspace##
bedtools intersect -v -a ${path1}/Extras/Mtruncatula_chromlen.bed -b ${BLACKLISTS}/Contigs_PlastidBL.bed >${GAT}/GAT_workspace.bed 

###Prepare Annotations###

#Homer: Not enough space annotated
#awk -F"\t" 'BEGIN {OFS="\t"} {print $2,$3,$4,$6,$6,$5}' ${ANNOTATION}/Mtruncatula_Homer_final.txt | sed 's/E/exon/g' | sed 's/I/intron/g' | sed 's/N/intergenic/g' | awk '$2<0 {$2=0} 1' OFS='\t' > ${ANNOTATION}/GAT1.bed
#bedtools intersect -f 0.5 -a ${ANNOTATION}/GAT1.bed -b ${GAT}/GAT_workspace.bed |sort -k1,1 -k2,2n  > ${GAT}/Homer.bed

echo "prepare annotation" 
#used:
bedtools intersect -f 0.5 -a ${ANNOTATION}/Mtruncatula_GAT.bed -b ${GAT}/GAT_workspace.bed > ${GAT}/Annot.bed
awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $8, $8, $6}' ${GAT}/Annot.bed > ${GAT}/GAT.bed

#loop (used for two Annotations, now obsolete#
Annot=(${GAT}/GAT.bed)

echo "prepare K27+ annotation"
bedtools intersect -f 0.5 -u -a  ${Annot} -b ${EPIC}/*intersect.bed | awk '{$4=$4 "_K27me3+"} {$5=$4} 1' OFS='\t'> ${GAT}/K27+_Annot_GAT.bed
bedtools intersect -f 0.5 -v -a ${Annot} -b  ${EPIC}/*intersect.bed | awk '{$4=$4 "_K27me3-"} {$5=$4} 1' OFS='\t'> ${GAT}/K27-_Annot_GAT.bed
cat ${GAT}/*Annot_GAT.bed | sort -k1,1 -k2,2n > ${GAT}/K27_all_GAT.bed

#Start faster/ other script
#bsub -q normal -R "rusage[mem=10240]" -M 12288 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/13_GAT*Kopie.sh

###Prepare workspaces###

#workspace for last 6 chromosomes
#tail -n 6 /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/GAT/GAT_workspace.bed > ${GAT}/GAT_workspace_6CHROMs.bed
#tail -n 12 ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_CHROMs.bed

#workspace minus Repeatregions + telobox rich regions (Blacklist) #
#gff2bed < ${DOWNLOADS}/*repeats*.gff > ${DOWNLOADS}/Aggressive_repeats.bed
#tail -n 12 ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_CHROMs.bed
#bedtools subtract -a GAT_workspace_CHROMs.bed -b ${DOWNLOADS}/Aggressive_repeats.bed > ${GAT}/GAT_workspace_norepeats.bed
#bedtools subtract -a ${GAT}/GAT_workspace_norepeats.bed -b ${BLACKLISTS}/Mtruncatula_telomerregions_BL.bed > ${GAT}/GAT_workspace_repeats+BL.bed
#bedtools subtract -a ${GAT}/GAT_workspace_CHROMs.bed -b ${ANNOTATION}/GFF3_intergenic.bed > ${GAT}/GAT_workspace_nointergenic.bed
#bedtools subtract -a ${GAT}/GAT_workspace_norepeats.bed -b ${ANNOTATION}/GFF3_intergenic.bed > ${GAT}/GAT_workspace_norpts+intergenic.bed
#bedtools subtract -a ${GAT}/GAT_workspace_repeats+BL.bed -b ${ANNOTATION}/GFF3_intergenic.bed > ${GAT}/GAT_workspace_norpts+intergenic+telorpts

#Loop for Gat analysis to use each chromosome by itself#
#for z in $(seq 1 $(wc -l < ${GAT}/GAT_workspace.bed))
#do
#sed -n ${z}p ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_CHROM.bed

#Loop for gat analysis to use two Chromosomes at once on workspace#
#for z in $(seq 1 $(wc -l < ${GAT}/GAT_workspace.bed))
#do
#sed -n ${z}p ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_1CHROM.bed
#sed -n $((${z}+1 | bc -l))p ${GAT}/GAT_workspace.bed > ${GAT}/GAT_workspace_2CHROM.bed
#cat ${GAT}/GAT_workspace_1CHROM.bed ${GAT}/GAT_workspace_2CHROM.bed > ${GAT}/GAT_workspace_2CHROMs.bed


### GAT Analyse ####

rm -f ${GAT}/gat.log

#divide Motif files by strain (only for strain dependent analysis)#
cd ${TELOBED}

#for i in *_Mtruncatula.bed
#do
#awk -F"\t" 'BEGIN {OFS="\t"} ($6=="+") {print $0}' $i > ${i%_Mtruncatula.bed}_+GATMtruncatula.bed
#awk -F"\t" 'BEGIN {OFS="\t"} ($6=="-") {print $0}' $i > ${i%_Mtruncatula.bed}_-GATMtruncatula.bed
#done/

#workspaces=("GAT_workspace_nointergenic.bed" ) #"GAT_workspace_norpts+intergenic.bed" "GAT_workspace_norpts+intergenic+telorpts"
#leftout=("intergenic" "rpts+intergenic" "telorpts+rpts+intergenic")

#for w in ${!workspaces[@]}
#do

echo "start GAT"

for i in *_Mtruncatula.bed
do
gat-run.py --ignore-segment-tracks \
	--segments=$i \
	--annotations=${GAT}/K27_all_GAT.bed \
	--workspace=${GAT}/GAT_workspace.bed\
	--num-samples=1000 \
	--log=${GAT}/gat.log > ${GAT}/${i%_Mtruncatula.bed}gat.out
#	--overlapping-annotations
awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%_Mtruncatula.bed}gat.out >${GAT}/${i%_Mtruncatula.bed}table.txt
awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%_Mtruncatula.bed}gat.out >${GAT}/${i%_Mtruncatula.bed}header.txt
cat ${GAT}/${i%_Mtruncatula.bed}header.txt ${GAT}/${i%_Mtruncatula.bed}table.txt > ${GAT}/${i%_Mtruncatula.bed}GAT_result.txt
done

#done

#done 
#for i in *_Mtruncatula.bed
#do
#gat-run.py --ignore-segment-tracks \
#	--segments=$i \
#	--annotations=${GAT}/K27_all_GAT.bed \
#	--workspace=${GAT}/GAT_workspace_2CHROMs.bed \
#	--num-samples=500 \
#	--log=${GAT}/gat.log > ${GAT}/${i%_Mtruncatula.bed}gat.out
#	--overlapping-annotations
#awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%_Mtruncatula.bed}gat.out >${GAT}/${i%_Mtruncatula.bed}table.txt
#awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%_Mtruncatula.bed}gat.out >${GAT}/${i%_Mtruncatula.bed}header.txt
#cat ${GAT}/${i%_Mtruncatula.bed}header.txt ${GAT}/${i%_Mtruncatula.bed}table.txt > ${GAT}/${i%_Mtruncatula.bed}GAT_result_CHR${z}_${z}+1.txt
#done
#done

#For strain dependent analysis #
#for i in *_-GATMtruncatula.bed
#do
#gat-run.py --ignore-segment-tracks \#
#	--segments=$i \
#	--annotations=${GAT}/GAT-.bed \
#	--workspace=${GAT}/GAT_workspace_CHROM.bed \
#	--num-samples=1000 \
#	--log=${GAT}/gat.log > ${GAT}/${i%GATMtruncatula.bed}gat.out
#	--overlapping-annotations
#awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%GATMtruncatula.bed}gat.out >${GAT}/${i%GATMtruncatula.bed}table.txt
#awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%GATMtruncatula.bed}gat.out >${GAT}/${i%GATMtruncatula.bed}header.txt
#cat ${GAT}/${i%GATMtruncatula.bed}header.txt ${GAT}/${i%GATMtruncatula.bed}table.txt > ${GAT}/${i%GATMtruncatula.bed}GAT_result_CHROM${z}.txt
#done
#done

#Run peaks vs Teloboxes: Failed#
#for i in *_Mtruncatula.bed
#do
#awk -F"\t" -vvar="${i%_Mtruncatula.bed}" 'BEGIN {OFS="\t"} {$4=var} {$5=$4}' $i > ${GAT}/${i%_Mtruncatula.bed}_Annot.bed
#done
#cat ${GAT}/*_Annot.bed > ${GAT}/Motif_GAT.bed
#gat-run.py --ignore-segment-tracks \
#	--segments=${EPIC}/*intersect_SigPeaks.bed \
#	--annotations=${GAT}/Motif_GAT.bed \
#	--workspace=${GAT}/GAT_workspace.bed \
#	--num-samples=100 \
#	--log=${GAT}/Motif_gat.log > ${GAT}/Motif_gat.out
#awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/Motif_gat.out >${GAT}/Motif_table.txt
#awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/Motif_gat.out >${GAT}/Motif_header.txt
#cat ${GAT}/Motif_header.txt ${GAT}/Motif_table.txt > ${GAT}/Motif_GAT_result.txt


rm -f ${GAT}/*table.txt
rm -f  ${GAT}/*header.txt
rm -f ${GAT}/*gat.out
rm -f ${GAT}/*1CHROM.bed ${GAT}/*2CHROM.bed ${GAT}/*2CHROMs.bed
