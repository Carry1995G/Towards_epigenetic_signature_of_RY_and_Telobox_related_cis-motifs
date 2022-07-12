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
