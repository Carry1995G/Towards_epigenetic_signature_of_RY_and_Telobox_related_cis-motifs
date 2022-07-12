#!/bin/bash

#requirement: python3.6
cd ~

GENES=/Chrom_FeatureOverlap
TELO_BIGWIG=/Motif_BIGWIG
MATRIX=/Matrix
BLACKLIST=/Blacklists
PLOT=/Plots

motifs=("Telobox" "Telolike" "RY")

##Matrix-Calculations ##

cd {TELO_BIGWIG}
for i in "${!motifs[@]}"
do
computeMatrix scale-regions -S ${motifs[i]}*.bw ${motifs[i]}*_shuffle* \
				-R ${GENES}/*+Genes.bed ${GENES}/*-Genes.bed \
				-b 1000 -a 1000 \
				-bl ${BLACKLIST}/Contigs*.bed \
				--missingDataAsZero \
				--skipZeros \
				-p 2 \
				-o ${MATRIX}/${motifs[i]}_H3K27me3.gz
done

##Generate a plot##
cd ${MATRIX}
for i in *_H3K27me3.gz
do
plotProfile -m $i \
		-o ${PLOT}/${i%.gz}.png \
		--samplesLabel "${i%_H3K27me3.gz} motifs" "shuffle control motifs" \
		--regionsLabel "H3K27me3+ genes" "H3K27me3- genes" \
		-y "Coverage per 50 bp" \
		--legendLocation "best" \
		--plotType se
done
