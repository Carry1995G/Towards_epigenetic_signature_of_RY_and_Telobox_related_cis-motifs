#!/bin/bash

#### Motif annotation using fuzznuc (EMBOSS) ####

cd /~

#New folders and path shortcuts
mkdir -p Motif_GFF
MOTIFGFF=/Motif_GFF

patterns=("AAACCCTA" "[GA][CA]CCTA[GA]" "CATGCA")
pnames=("Telobox" "Telolike" "RY")
for p in "${!patterns[@]}"
do
fuzznuc ./GENOME.fasta \
	-pattern "${patterns[p]}" \
	-complement \
	-rformat gff ./Motif_GFF/OUTPUT.gff
done

#Converting annotated GFF-files to BED and BAM-files#

#New folders and path shortcuts
mkdir -p Motif_BAM
mkdir -p Motif_BED
mkdir -p Extras
EXTRAS=/Extras
BED=/Motif_BED
BAM=/Motif_BAM
BLACKLIST=/Blacklists #if available

#Generate files containing chromosome length (different formats)#
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

#Generate shuffled motif BED- and BAM-files for negative control#
cd ${BED}
for i in OUTPUT.bed
do
bedtools shuffle -i $i \
	-g ${EXTRAS}/Chromlen.chromsizes > ${BED}/${i%OUTPUT.bed}shuffle.bed	
bedtools bedtobam -i ${BED}/${i%OUTPUT.bed}shuffle.bed \
	-g ${EXTRAS}/Chromlen.chromsizes | samtools sort -o ${BAM}/${i%OUTPUT.bed}shuffle.bam -
samtools index ${BAM}/${i%OUTPUT.bed}shuffle.bam
done

#### Motif coverage calculation ####

#New folders and path shortcuts
mkdir -p Motif_BIGWIG
MOTIFBW=/Motif_BIGWIG

cd ${BAM}
for i in *.bam
do
bamCoverage -b $i -o $[MOTIFBW}/${i%.bam}.bw -p 70 --extendReads 200 --ignoreDuplicates -bl ${BLACKLIST}/Contigs*.bed
done
