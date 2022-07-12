#!/bin/bash

### Motif Annotation mithilfe von fuzznuc (EMBOSS)
#bsub -q normal -R 10240 -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/3_Motif_annotation.sh

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/
mkdir -p Motif_GFF

patterns=("AAACCCTA" "[GA][CA]CCTA[GA]" "CATGCA")
pnames=("Telobox" "Telolike" "RY")

for p in "${!patterns[@]}"
do
fuzznuc ./Downloads/*genome.fasta \
	-pattern "${patterns[p]}" \
	-complement \
	-rformat gff ./Motif_GFF/${pnames[p]}_Mtruncatula.gff
done

bsub -q bigmem -R "rusage[mem=20000]" -M 40000 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/4_*