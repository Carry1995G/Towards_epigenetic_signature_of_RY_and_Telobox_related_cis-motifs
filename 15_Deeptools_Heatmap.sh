#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/15_Deeptools_He*

###conda activate py36
cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

GENES=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Chrom_FeatureOverlap
TELO_BIGWIG=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Motif_BIGWIG
MATRIX=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Matrix
BLACKLIST=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists
PLOT=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Plots


motifs=("Telobox" "Telolike" "RY")

###Matrix-Calculations ###
echo "Calculate Matrix"

#computeMatrix scale-regions -S ${TELO_BIGWIG}/Telobox*.bw \
#				${TELO_BIGWIG}/Telolike*.bw \
#				-R ${GENES}/*merged*K27+Genes.bed \
#				${GENES}/*merged*K27-Genes.bed \
#				-b 1000 -a 1000 \
#				--missingDataAsZero \
#				-bl ${BLACKLIST}/Min_Blacklist.bed \
##				--skipZeros \
#				-p 2 \
#				-o ${MATRIX}/Telobox_Telolike_minBL_K27+-.gz

###Heatmap-Calculations ###

#plotHeatmap 	-m ${MATRIX}/Telobox_Telolike_minBL_K27+-.gz \
#		-o ${PLOT}/Telobox_Telolike_minBL_27+-.png \
#		--samplesLabel "Telobox" "Telolike" \
#		--regionsLabel "H3K27+ Genes" "H3K27- Genes" \
#		-y "Number of motifs per 50 bp" \
#		--legendLocation "best"

#########used:
cd ${TELO_BIGWIG}

###Matrix-Calculations ###
echo "Calculate Matrix"

##SYNTAX:##
#computeMatrix scale-regions -S MOTIF-BIGWIG \
#				-R GENES +/- H3K27me3 \
#				-b 1000 -a 1000 \
#				--missingDataAsZero \
#				-bl BLACKLIST \
##				--skipZeros \
#				-p 2 \
#				-o GZ-FILE


for i in "${!motifs[@]}"
do
computeMatrix scale-regions -S ${motifs[i]}*_Mtruncatula.bw ${motifs[i]}*_shuffleM* \
				-R ${GENES}/*+Genes.bed ${GENES}/*-Genes.bed \
				-b 1000 -a 1000 \
				-bl ${BLACKLIST}/Contigs*.bed \
				--missingDataAsZero \
				--skipZeros \
				-p 2 \
				-o ${MATRIX}/${motifs[i]}_Mtruncatula_BL15kb_K27+-shuffle.gz
done

#for i in "${!motifs[@]}"
#do
#computeMatrix scale-regions -S ${motifs[i]}*_Mtruncatula.bw ${motifs[i]}*_shuffleM* \
#				-R ${GENES}/*+CDS.bed ${GENES}/*-CDS.bed \
#				-b 1500 -a 1500 \
#				-bl ${BLACKLIST}/Contigs*.bed \
#				--missingDataAsZero \
#				--skipZeros \
#				-p 2 \
#				-o ${MATRIX}/${motifs[i]}_Mtruncatula_BL3kb_CDSK27+-shuffle.gz
#done


#for i in *_Mtruncatula.bw
#do
#computeMatrix scale-regions -S $i \
#				-R ${GENES}/*+Genes.bed ${GENES}/*-Genes.bed \
#				-b 1000 -a 1000 \
#				--missingDataAsZero \
#				--skipZeros \
#				-p 2 \
#				-o ${MATRIX}/${i%.bw}_K27+-.gz
#done

###Heatmap-Calculations ###
echo "Generate Heatmap"

cd ${MATRIX}

for i in *Mtruncatula_BL15kb_K27+-shuffle.gz
do
plotProfile -m $i \
		-o ${PLOT}/${i%.gz}.png \
		--samplesLabel "${i%_Mtruncatula_BL15kb_K27+-shuffle.gz} motifs" "shuffle control motifs" \
		--regionsLabel "H3K27me3+ genes" "H3K27me3- genes" \
		-y "Coverage per 50 bp" \
		--legendLocation "best" \
		--plotType se
done

#for i in *_Mtruncatula_BL3kb_CDSK27+-shuffle.gz
#do
#plotProfile -m $i \
#		-o ${PLOT}/${i%.gz}.png \
#		--samplesLabel "${i%_Mtruncatula_BL_CDSK27+-shuffle.gz}-Motif" "shuffle control motif" \
#		--regionsLabel CDS-H3K27me3+ CDS-H3K27me3- \
#		-y "Number of motifs per 50 bp" \
#		--legendLocation "best" \
#		--plotType se
#done


#for i in *Mtruncatula_K27+-.gz
#do
#plotProfile -m $i \#
#		-o ${PLOT}/${i%.gz}_$p.png \
#		--samplesLabel "${i%_Mtruncatula_K27+-.gz}-Motif" \
#		--regionsLabel H3K27me3+ H3K27me3- \
#		-y "Number of motifs per 50 bp" \
#		--legendLocation "best" \
#		--plotType se
#done