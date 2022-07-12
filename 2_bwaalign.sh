#!/bin/bash

#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/2_bwa*
#bsub -q bigmem -R "rusage[mem=20000]" -M 40000 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/2_bwa*

#######Create bwa index for tomato. Needs to be done only once############

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/
mkdir -p Downloads
mkdir -p BWA-align
mkdir -p BWA-align/BAM-files
mkdir -p BWA-align/BED-files

#Paths
BWA=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/BWA-align/
DOWNLOADS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Downloads

cd ${DOWNLOADS}
#not used for the SICER-sequences 
#MtrunA17r5.0-ANR genome files GCA
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/473/485/GCA_003473485.2_MtrunA17r5.0-ANR/GCA_003473485.2_MtrunA17r5.0-ANR_genomic.fna.gz
#gzip -d GCA_003473485.2_MtrunA17r5.0-ANR_genomic.fna.gz

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/473/485/GCA_003473485.2_MtrunA17r5.0-ANR/GCA_003473485.2_MtrunA17r5.0-ANR_genomic.gff.gz
#gzip -d GCA_003473485.2_MtrunA17r5.0-ANR_genomic.gff.gz

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/473/485/GCA_003473485.2_MtrunA17r5.0-ANR/GCA_003473485.2_MtrunA17r5.0-ANR_genomic.gtf.gz
#gzip -d GCA_003473485.2_MtrunA17r5.0-ANR_genomic.gtf.gz

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/473/485/GCF_003473485.1_MtrunA17r5.0-ANR/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna.gz
#gzip -d GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna.gz

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/473/485/GCF_003473485.1_MtrunA17r5.0-ANR/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gtf.gz
#gzip -d GCF_003473485.1_MtrunA17r5.0-ANR_genomic.gtf.gz
#MtrunA17r5.0-ANR genome files GCF

#used:
wget https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.8/MtrunA17r5.0-ANR-EGN-r1.8.fastaFiles.zip
unzip MtrunA17r5.0-ANR-EGN-r1.8.fastaFiles.zip

wget https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.8/MtrunA17r5.0-ANR-EGN-r1.8.gtf.zip
unzip MtrunA17r5.0-ANR-EGN-r1.8.gtf.zip

wget https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.8/MtrunA17r5.0-ANR-EGN-r1.8.structuralAnnotation.zip
unzip MtrunA17r5.0-ANR-EGN-r1.8.structuralAnnotation.zip

#For GO-Analysis
wget https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.8/MtrunA17r5.0-ANR-EGN-r1.8.functionalAnnotation.20210505.zip
unzip MtrunA17r5.0-ANR-EGN-r1.8.functionalAnnotation.20210505.zip

#gzip Pecrix*.zip

##Not used due to SICER H3K27me3 files!!

#bwa index -p ${BWA}/Mtruncatulabwaidx -a bwtsw GCF_000004515.6_Glycine_max_v4.0_genomic.fna

#cd ../fastq

#for i in "SRR9622334" "SRR9622335" "SRR9622336" "SRR9622337" "SRR9622338"
#do
#bwa mem -t 4 ${BWA}/Mtruncatulabwaidx ${i}_1.fastq ${i}_2.fastq | \
#			samtools view -b - | \
#			samtools sort - -o ${BWA}/BAM-files/${i}.bam
#done


####make index files for BAM files####
#cd ${BWA}/BAM-files

#mv SRR9622334.bam Input-R1.bam 
#mv SRR9622335.bam Input-R2.bam 
#mv SRR9622336.bam H3K27me3-R1.bam 
#mv SRR9622337.bam H3K27me3-R2.bam 
#mv SRR9622338.bam H3K27me3-R3.bam 

####Shuffle control files ####
#holt 50% der Reads raus (2 ist der random generator .50 die Fraktion)
#-h muss dabei sein, damit der header beleibt.
#cd ${BWA}/BAM-files
#samtools view -b -h -s 2.50 -b Input-R1.bam > Input-R1_50p.bam
#samtools view -b -h -s 2.50 -b Input-R2.bam > Input-R2_50p.bam
#samtools merge Input-merged.bam Input-R1_50p.bam Input-R2_50p.bam 
#rm Input*50p.bam

#for i in *.bam
#do
#samtools index $i >${i%.bam}.bai
#done

#####make bedfiles from the bam files######
#for i in *.bam
#do
#bedtools bamtobed -i $i > ${BWA}/BED-files/"${i%.bam}.bed"
#done


#####not done####
#count the mapped reads for each mapping

#cd ${BWA}/BED-files
#for i in *bed
#do
#wc -l $i >>mapped_reads.txt
#done

#bsub -q normal -R 10240 -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/3_Motif_annotation.sh
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/9_EPIC_Peakcalling.sh