#!/bin/bash

### Motif annotation using fuzznuc (EMBOSS)###
cd /~
mkdir -p Motif_GFF

patterns=("AAACCCTA" "[GA][CA]CCTA[GA]" "CATGCA")
pnames=("Telobox" "Telolike" "RY")
for p in "${!patterns[@]}"
do
fuzznuc ./INPUT.fasta \
	-pattern "${patterns[p]}" \
	-complement \
	-rformat gff ./Motif_GFF/OUTPUT.gff
done


### Motif coverage calculation ###

cd /~

mkdir -p Motif_BAM
mkdir -p Motif_BED
mkdir -p Extras
mkdir -p Motif_BIGWIG

#Path shortcuts
DOWNLOADS=/Downloads
EXTRAS=/Extras
BED=/Motif_BED
BAM=/Motif_BAM
BLACKLIST=/Blacklists

#Generate files containing chromosome length (different formats)#
cd ${Extras}
pip install pyfaidx
faidx ${DOWNLOADS}/GENOME-FASTAFILE.fasta -i chromsizes > Chromlen.txt #used for later steps
faidx ${DOWNLOADS}/GENOME-FASTAFILE.fasta -i chromsizes > Chromlen.bed #used for later steps
cp Chromlen.txt Chromlen.chromsizes

#Reformat annotated GFF-files to BAM-files#
cd ../Motif_GFF
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

#Coverage calculation of all regions #
cd ${BAM}
for i in *.bam
do
bamCoverage -b $i -o ../Motif_BIGWIG/${i%.bam}.bw -p 70 --extendReads 200 --ignoreDuplicates -bl ${BLACKLIST}/Contigs*.bed
done
