#!/bin/bash

#### Motif and negativ control annotation####

##Motif annotation##
cd /~

#new folders and path shortcuts
mkdir -p Motif_GFF
MOTIFGFF=/Motif_GFF

#annotate all three motifs in a single loop using fuzznuc (EMBOSS) and save in three seperate gff files#
#The motif names are saved in the beginning of the file title and remain there throughout the entire analysis#

patterns=("AAACCCTA" "[GA][CA]CCTA[GA]" "CATGCA")
pnames=("Telobox" "Telolike" "RY")
for p in "${!patterns[@]}"
do
fuzznuc ./GENOME.fasta \
	-pattern "${patterns[p]}" \
	-complement \
	-rformat gff ./Motif_GFF/${pnames[p]}_OUTPUT.gff
done

##Convert annotated GFF-files to BED and BAM-files##

#New folders and path shortcuts
mkdir -p Motif_BAM
mkdir -p Motif_BED
mkdir -p Extras
EXTRAS=/Extras
BED=/Motif_BED
BAM=/Motif_BAM
BLACKLIST=/Blacklists #if available

#Generate files containing chromosome length in different file formats using pyfaidx#
cd ${Extras}
pip install pyfaidx
faidx GENOME.fasta -i chromsizes > Chromlen.txt #used for later steps
faidx GENOME.fasta -i chromsizes > Chromlen.bed #used for later steps
cp Chromlen.txt Chromlen.chromsizes

#Reformat annotated GFF-files to BAM-files#
cd ${MotifGFF}
for i in *OUTPUT.gff
do
sed -i '/^[#]/ d' $i 			#remove fuzznuc comment 
sed -i '/^$/d' $i 			#remove any space 
gff2bed < $i | bedtools bedtobam -i - -g ${EXTRAS}/Chromlen.chromsizes | samtools sort -o ${BAM}/"${i%.gff}.bam" -
samtools index ${BAM}/${i%.gff}OUTPUT.bam
gff2bed < $i > ${BED}/"${i%.gff}OUTPUT.bed"
done

##Obtain negative control by shuffling motif positions of BED file and save as BED and BAM-files##
cd ${BED}
for i in *OUTPUT.bed
do
bedtools shuffle -i $i \
	-g ${EXTRAS}/Chromlen.chromsizes > ${BED}/${i%OUTPUT.bed}shuffle.bed	
bedtools bedtobam -i ${BED}/${i%OUTPUT.bed}shuffle.bed \
	-g ${EXTRAS}/Chromlen.chromsizes | samtools sort -o ${BAM}/${i%OUTPUT.bed}shuffle.bam -
samtools index ${BAM}/${i%OUTPUT.bed}shuffle.bam
done

#### Motif coverage calculation and visualization ####

#Motif coverage calculation with deeptools bamCoverage#

#New folders and path shortcuts
mkdir -p Motif_BIGWIG
MOTIFBW=/Motif_BIGWIG

cd ${BAM}
for i in *.bam
do
bamCoverage -b $i -o $[MOTIFBW}/${i%.bam}.bw -p 70 --extendReads 200 --ignoreDuplicates -bl ${BLACKLIST}/Contigs*.bed
done

####Visualization####

cd ~
#new folders and path shortcut
mkdir -p Matrix
mkdir -p Plots
PLOTS=/Plots
MATRIX=/Matrix

##Compute Matrix with deeptools computeMatrix ##
cd ${MOTIFBW}

#Generate matrix for each motif (except shuffled motifs as negative control)#
pnames=("Telobox" "Telolike" "RY")
for i in *OUTPUT.bw						#only the files of the real motif positioned are selected, not the shuffled positions
do
computeMatrix scale-regions \
				-S $i \
				-R ./ANNOTATION.bed \
				-b 1000 -a 1000 \
				-bl ${BLACKLIST}/Contigs*.bed \
				-o ${MATRIX}/${i%.bw}MATRIX.gz
done

##Generate plot profiles of the matrices with deeptools plotProfile## 
cd Matrix

for i in *MATRIX.gz
do
plotProfile -m $i  \
			-o ${PLOTS}/${i%MATRIX.gz}PLOT.png \
			--samplesLabel "${i%MATRIX.gz} motifs" \ #removes end of file name so that only motif name remains
			--plotType se \
			-y "Coverage per 50 bp" \
			--legendLocation "best"
done
