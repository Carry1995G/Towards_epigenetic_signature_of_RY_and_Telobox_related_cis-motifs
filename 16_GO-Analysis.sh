#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/16_GO*

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

#https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/web/fatal/entry.html?accession=MtrunA17_Chr1g0147081

mkdir -p GO-Analysis

GENES=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Chrom_FeatureOverlap
TELO_BIGWIG=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Motif_BIGWIG
MATRIX=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Matrix
BLACKLIST=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists
PLOT=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Plots
GO=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/GO-Analysis
ANNOTATION=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Annotation_Files


motifs=("Telobox" "Telolike" "RY")
cd ${TELO_BIGWIG}

for i in "${!motifs[@]}"
do
computeMatrix scale-regions -S ${motifs[i]}*_Mtruncatula.bw \
				-R ${GENES}/*+Genes.bed \
				-bl ${BLACKLIST}/Contigs_PlastidBL.bed \
				-b 1000 -a 1000 \
				--missingDataAsZero \
				--skipZeros \
				-p 2 \
				-o ${MATRIX}/${motifs[i]}_Mtruncatula_GenesGOkmeans.gz
done

cd ${MATRIX}

for i in "${!motifs[@]}"
do
plotProfile -m ${motifs[i]}_Mtruncatula_GenesGOkmeans.gz \
		-y "Coverage per 50 bp" \
		--samplesLabel "${motifs[i]} motifs" \
		--legendLocation "best" \
		--kmeans 5\
		--plotType heatmap \
		-o ${PLOT}/${motifs[i]}_Mtruncatula_GenesGOkmeans.png \
	--outFileSortedRegions ${GO}/${motifs[i]}_Mtruncatula_GenesGOkmeans.txt
done


cd ${GO}

#how the files look like:
#MtrunA17Chr1	1417972	1418204	Genes:MtrunA17Chr1g0147901.1	.	+	1417972	1418204	0	1	232	1417971	cluster_1
#split $4 at ":" and retrieve second party with awk, then insert _ after MtrunA17 with sed (GO-files use this annotation)
#for i in *GenesGOkmeans.txt
#do 
#awk -vFS="\t" -vOFS="\t" ' split($4, a, /[:]/) >= 2 {print a[2], $13}' $i | sed 's/MtrunA17/&_/' - > ${i%_GenesGOkmeans.txt}_GenesinClusters.txt
#done

for i in *_GenesGOkmeans.txt
do 
awk -vFS="\t" -vOFS="\t" ' split($4, a, /[:]/) >= 2 {print a[2], $13}' $i | sed 's/MtrunA17/&_/' - > ${i%_GOkmeans.txt}_GenesinClusters.txt
done

#for x in *_GenesinClusters.txt
#do
#sed 's/MtrunA17/&_/4' $x
#done
#awk -vFS="\t" -vOFS="\t" '$3 ~ /cluster_1|cluster_2/ {print $2}' Telobox_GenesinClusters.txt >Telobox_Genes_enriched_TSS.txtc


#set up list of genes 
awk -vFS="\t" -vOFS="\t" 'split($4, a, /[:]/) >= 2 {print a[2]}' ${ANNOTATION}/GFF2BED_genes.bed | sed 's/MtrunA17/&_/' - | awk /MtrunA17_Chr[1-9]/ > Mtruncatula_Genes.txt