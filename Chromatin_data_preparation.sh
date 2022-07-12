#!/bin/bash

### Alignment of chromatin data to Genome###
#new folders and paths
mkdir -p BWA-align
mkdir -p BWA-align/BAM-files
mkdir -p BWA-align/BED-files
BWA=/BWA-align
DOWNLOADS=/Downloads

#Generate bwa index#
bwa index -p ${BWA}/bwaidx -a bwtsw GENOMIC.fasta

#Align to genome#
#cd ../fastq

Depending on the sequencing method used on the downloaded data (paired-end or single read), a different loop needs to be chosen#

#Paired end data:#
SRRfiles=("SRRxxx" "SRRyyy" "SRRzzz")
for i in ${SRRfiles}
do
bwa mem -t 4 ${BWA}/waidx ${i}_1.fastq ${i}_2.fastq | \
			samtools view -b - | \
			samtools sort - -o ${BWA}/BAM-files/${i}.bam
done

#Single-read data:# 
for i in *.fastq
do
bwa mem -t 4 ${BWA}/bwaidx ${i} | \
			samtools view -b - | \
			samtools sort - -o ${BWA}/BAM-files/${i}.bam
done

#make index files for BAM files and convert to BED files#

cd ${BWA}/BAM-files
for i in *.bam
do
samtools index $i >${i%.bam}.bai \
bedtools bamtobed -i $i > ${BWA}/BED-files/"${i%.bam}.bed"
done


###Identification of enriched genes and regions###

mkdir -p Epic
EPIC=/Epic

##Identification of H3K27me3 regions: Peak calling with EPIC##


#Save overlapping H3K27me3 positions of Replicates#
bedtools intersect -a ${H3K27ME3}1/*R1 -b ${H3K27ME3}2/*R2 |  \ 			#intersection of two EPIC files
bedtools intersect -wa -a - -b ${EXTRAS}/CHROMLEN.bed | \				#intersection with Genome to remove erronous positions
sort -k1,1 -k2,2n > ${PEAKS}/H3K27me3_intersect.bed					#sort

##Identification of H3K27me3 positive and negative genes##

#Filter gene positions:
gff2bed < ${Downloads}/GENOME.GFF| awk '$8=="gene" {print $0}'> ${Annot}/GFF2BED_genes.bed

cd ${EPIC}

#Intersect H3K27me3-positions with gene positions and save genes overlapping with H3K27me3
for i in *intersect.bed
do
bedtools intersect -a ${Annot}/GFF2BED_genes.bed -b $i -f 0.5 -u > ${Genes}/${i%me3_intersect.bed}K27+Genes.bed \
awk -vFS="\t" -vOFS="\t" 'split($4, a, /[:]/) >= 2 {print a[2]}' ${Genes}/${i%me3_intersect.bed}K27+Genes.bed | \
sed 's/MtrunA17/&_/' - | awk /MtrunA17_Chr[1-9]/ | uniq > ${GENE}/${i%me3_intersect.bed}K27+Genes.txt
done

#Intersect H3K27me3-positions with gene positions and save genes not overlapping with Hh3K27me3
for i in *intersect.bed
do
bedtools intersect -a ${Annot}/GFF2BED_genes.bed -b $i -f 0.5 -v > ${Genes}/${i%me3_intersect.bed}K27-Genes.bed \
done
