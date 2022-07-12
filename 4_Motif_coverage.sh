#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/4_*
#bsub -q bigmem -R "rusage[mem=20000]" -M 40000 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/4_*

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

mkdir -p Motif_BAM
mkdir -p Motif_BED
mkdir -p Extras
mkdir -p Motif_BIGWIG

#BLACKLISTS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists
DOWNLOADS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Downloads
EXTRAS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Extras
BED=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Motif_BED
BAM=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Motif_BAM
BLACKLIST=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists


####Make fasta chromosome length file ###
cd Extras
echo "Chromosomlängen-Datei wird hergestellt..."

#awk '$0 ~ ">" OFS='\t' {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' ${DOWNLOADS}/S_lycopersicum*.fa > Tomato_chromlen.txt
#sed -i 's/Solanum lycopersicum cultivar Heinz 1706 unplaced genomic scaffold, SL3.0 SL3.00SC000//g' Tomato_chromlen.txt
#sed -i 's/Solanum lycopersicum cultivar Heinz 1706 chromosome [0-9], SL3.0, whole genome shotgun seque//g' Tomato_chromlen.txt
#sed -i 's/Solanum lycopersicum cultivar Heinz 1706 chromosome [0-9]*, SL3.0, whole genome shotgun sequ//g' Tomato_chromlen.txt
#sed -i 's/Solanumlycopersicumbio-materialTGRC:LA1421mitochondrion,completegenome//g' Tomato_chromlen.txt
#sed -i 's/Solanumlycopersicumchloroplast,completegenome//g' Tomato_chromlen.txt
#sed -i 's/ //g' Tomato_chromlen.txt

#pip install pyfaidx
#faidx ${DOWNLOADS}/*genome.fasta -i chromsizes > Mtruncatula_chromlen.txt
#faidx ${DOWNLOADS}/*genome.fasta -i chromsizes > Mtruncatula_chromlen.bed

cp Mtruncatula_chromlen.txt Mtruncatula.chromsizes


###Reformat annotated GFF-files to BAM-files##
#echo "Umformatieren der annotierten GFF-Files in BED und BAM Files..."

cd ../Motif_GFF

#for i in *Mtruncatula.gff
#do
#sed -i '/^[#]/ d' $i #damit wird man die Kommentarzeile des outputs von fuzznuc los
#sed -i '/^$/d' $i # damit wird man Leerzeilen los
#gff2bed < $i | bedtools bedtobam -i - -g ${EXTRAS}/Mtruncatula.chromsizes | samtools sort -o ../Motif_BAM/"${i%.gff}.bam" -
#samtools index ${BAM}/${i%.gff}.bam
#gff2bed < $i > ${BED}/"${i%.gff}.bed"
#done

#Shuffle BED-files for Control
cd ${BED}

#for i in *Mtruncatula.bed
#do
#bedtools shuffle -i $i \
#	-g ${EXTRAS}/Mtruncatula.chromsizes > ${BED}/${i%Mtruncatula.bed}shuffleMtruncatula.bed #only chromosome length needed (ch1 xxxx), not chromosomone size (ch1 0 xxxx) 
#bedtools bedtobam -i ${BED}/${i%Mtruncatula.bed}shuffleMtruncatula.bed \
#	-g ${EXTRAS}/Mtruncatula.chromsizes| samtools sort -o ${BAM}/${i%Mtruncatula.bed}shuffleMtruncatula.bam -
#samtools index ${BAM}/${i%Mtruncatula.bed}shuffleMtruncatula.bam
#done

###Coverage calculation of all regions ###
echo "Coverage wird berechnet..:"

cd ${BAM}

for i in *Mtruncatula.bam
do
bamCoverage -b $i -o ../Motif_BIGWIG/${i%.bam}.bw -p 70 --extendReads 200 --ignoreDuplicates -bl ${BLACKLIST}/Contigs*.bed
done

#bsub -q normal -R "rusage[mem=10240]" -M 12288 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/13_GAT.sh
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/6_AnnotGTF*
bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/7_ComputeMatrix*
bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/15_Deeptools_He*