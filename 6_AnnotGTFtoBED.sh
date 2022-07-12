#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/6_AnnotGTF*
###Change GTF into 6 column BED format 
cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

mkdir -p /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Annotation_Files
ANNOT=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Annotation_Files
cd ${ANNOT} 

#extract gene names and positions from annotation


awk '$3=="gene" {print $1"\t"$4-1"\t"$5"\t"$9"\t"$6"\t"$7}' ../Downloads/Mtrun*.gff3 | \
awk -F'\t' -v OFS='\t' ' split($4, a, /[:;]/) >= 2 {print $1, $2, $3, a[2], $5, $6}' > ${ANNOT}/Mtruncatula_6col.bed 
awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$9,$6,$7}}' ../Downloads/Mtrun*.gff3 | tr -d '";' >Mtruncatula.bed


bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/7_ComputeMatrix*