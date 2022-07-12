#!/bin/bash
##Compute Matrix for later profile plots ##

#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/7_ComputeMatrix*

BLACKLIST=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists


cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/
mkdir -p Matrix
cd Motif_BIGWIG

#used, but not needed again: 
#sed -i '/^[#]/ d' ../Downloads/*.gtf  #Kommentarzeile des outputs von fuzznuc loswerden

echo "Compute Matrix "

pnames=("Telobox" "Telolike" "RY")

###Minimale Blacklist
for i in *_Mtruncatula.bw
do
computeMatrix scale-regions \
				-S $i \
				-R ../Annotation_Files/*6col.bed \
				-b 1000 -a 1000 \
				-bl ${BLACKLIST}/Contigs*.bed \
				-o /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Matrix/${i%bw}gz
done

#computeMatrix scale-regions \
#				-S RY_Mtruncatula.bw Telobox_Mtruncatula.bw Telolike_Mtruncatula.bw\
#				-R ../Annotation_Files/*6col.bed \
#				-b 1000 -a 1000 \
#				-o /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Matrix/Motifs_Mtruncatula.gz

#Full Blacklist 

#for i in Telo*Tomato.bw
#do
#computeMatrix scale-regions \
#				-S $i \
#				-R ../Annotation_Files/*6col.bed \
#				-bl ../Blacklists/Telomer* \
#				-b 1000 -a 1000 \
#				-o /netscratch/dep_coupland/grp_turck/carina/Tomato/Matrix/${i%bw}_TeloChrBL.gz
#done

bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/8_Deeptools_Plot*