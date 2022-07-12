#!/bin/bash

#variables and paths
cd ~
mkdir -p GO-Analysis
GENES=/Chrom_FeatureOverlap
TELO_BIGWIG=/Motif_BIGWIG
MATRIX=/Matrix
BLACKLIST=/Blacklists
PLOT=/Plots
GO=/GO-Analysis
ANNOTATION=/Annotation_Files


motifs=("Telobox" "Telolike" "RY")

##Generate Matrix of Motifs##
cd ${TELO_BIGWIG}
for i in "${!motifs[@]}"
do
computeMatrix scale-regions -S ${motifs[i]}*.bw \
				-R ${GENES}/*+Genes.bed \
				-bl ${BLACKLIST}/Contigs_PlastidBL.bed \
				-b 1000 -a 1000 \
				--missingDataAsZero \
				--skipZeros \
				-p 2 \
				-o ${MATRIX}/${motifs[i]}GOkmeans.gz
done

##Genearte Plot from Matrix, thereby cluster with kmeans##
cd ${MATRIX}
for i in "${!motifs[@]}"
do
plotProfile -m ${motifs[i]}GOkmeans.gz \
		-y "Coverage per 50 bp" \
		--samplesLabel "${motifs[i]} motifs" \
		--legendLocation "best" \
		--kmeans 5\
		--plotType heatmap \
		-o ${PLOT}/${motifs[i]}GOkmeans.png \
	--outFileSortedRegions ${GO}/${motifs[i]}GOkmeans.txt
done


#print Gene name and cluster for each gene
cd ${GO}

for i in *_GenesGOkmeans.txt
do 
awk -vFS="\t" -vOFS="\t" '{print $4, $13}' $i > ${i%_GOkmeans.txt}_GenesinClusters.txt
done

#set up list of genes 
awk -vFS="\t" -vOFS="\t" 'split($4, a, /[:]/) >= 2 {print a[2]}' ${ANNOTATION}/GFF2BED_genes.bed | sed 's/MtrunA17/&_/' - | awk /MtrunA17_Chr[1-9]/ > Mtruncatula_Genes.txt
