#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Tomato/13_GAT*Kopie.sh
#bsub -q bigmem -R "rusage[mem=20000]" -M 40000 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Tomato/13_GAT*Kopie.sh

cd /netscratch/dep_coupland/grp_turck/carina/Tomato/

####################################################################
#variables and paths
path1=$(pwd)
mkdir -p ${path1}/GAT
GAT=${path1}/GAT/GAT_Fast
Gat=${path1}/GAT
IDR=${path1}/IDR
TELOBED=${path1}/Motif_BED
ANNOTATION=/netscratch/dep_coupland/grp_turck/carina/Tomato/Annotation_Files
Downloads=/netscratch/dep_coupland/grp_turck/carina/Tomato/Downloads


##Workspace##
bedtools subtract -a ${path1}/Extras/Tomato_chromlen.bed -b ${path1}/Blacklists/Min_Blacklist.bed > /netscratch/dep_coupland/grp_turck/carina/Tomato/GAT/GAT_workspace.bed

##Annotations##
#irgendwas doppelt sich hier (evtl gelöst)
#OLD

#awk -F"\t" 'BEGIN {OFS="\t"} {print $2,$3,$4,$6,$6,$5}' ${ANNOTATION}/Tomato_Homer_final.custom.txt | sed 's/E/exon/g' | sed 's/I/intron/g' | awk '$2<0 {$2=0} 1' OFS='\t' > ${ANNOTATION}/GAT1.bed
#sort -k1,1 -k2,2n ${ANNOTATION}/GAT1.bed > ${ANNOTATION}/sorted_GAT1.bed
#bedtools intersect -f 0.5 -a ${ANNOTATION}/sorted_GAT1.bed -b ${GAT}/GAT_workspace.bed  > ${GAT}/GAT.bed
#bedtools intersect -f 0.5 -u -a ${GAT}/GAT.bed -b ${IDR}/*merged*SigpValue.bed | awk '{$4=$4 "_K27+"} {$5=$4} 1' OFS='\t'> ${GAT}/K27+_Annot_GAT.bed
#bedtools intersect -f 0.5 -v -a ${GAT}/GAT.bed -b  ${IDR}/*merged*SigpValue.bed | awk '{$4=$4 "_K27-"} {$5=$4} 1' OFS='\t'> ${GAT}/K27-_Annot_GAT.bed
#cat ${GAT}/*Annot_GAT.bed > ${GAT}/K27_all_GAT.bed


#bedtools intersect -f 0.5 -a ${ANNOTATION}/Tomato_GAT.bed -b ${GAT}/GAT_workspace.bed  > ${GAT}/Annot.bed
#awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $8, $8, $6}' ${GAT}/Annot.bed | sort -k1,1 -k2,2n > ${GAT}/GAT.bed
#bedtools intersect -f 0.5 -u -a ${GAT}/GAT.bed -b ${IDR}/*merged*SigpValue.bed | awk '{$4=$4 "_K27+"} {$5=$4} 1' OFS='\t'> ${GAT}/K27+_Annot_GAT.bed
#bedtools intersect -f 0.5 -v -a ${GAT}/GAT.bed -b  ${IDR}/*merged*SigpValue.bed | awk '{$4=$4 "_K27-"} {$5=$4} 1' OFS='\t'> ${GAT}/K27-_Annot_GAT.bed
#cat ${GAT}/*Annot_GAT.bed | sort -k1,1 -k2,2n > ${GAT}/K27_all_GAT.bed


#######################

sed -n 1p /netscratch/dep_coupland/grp_turck/carina/Tomato/GAT/GAT_workspace.bed > /netscratch/dep_coupland/grp_turck/carina/Tomato/GAT/GAT_workspace_contigs.bed
tail -n 12p /netscratch/dep_coupland/grp_turck/carina/Tomato/GAT/GAT_workspace.bed > /netscratch/dep_coupland/grp_turck/carina/Tomato/GAT/GAT_workspace_CHROMs.bed

rm ${GAT}/gat.log

bedtools subtract -a ${Gat}/GAT_workspace_norepeats.bed -b ${BLACKLISTS}/Tomato_telomerregions_BL.bed > ${Gat}/GAT_workspace_repeats+BL.bed


cd ${TELOBED}
for i in *_Tomato.bed
do
for c in "contigs" "CHROMs"
do
gat-run.py --ignore-segment-tracks \
	--segments=$i \
	--annotations=/netscratch/dep_coupland/grp_turck/carina/Tomato/GAT/K27_all_GAT.bed \
	--workspace=/netscratch/dep_coupland/grp_turck/carina/Tomato/GAT/GAT_workspace_repeats+BL.bed \
	--num-samples=1000 \
	--log=${GAT}/gat.log > ${GAT}/${i%Tomato.bed}gat.out
awk -F"\t" 'BEGIN {OFS="\t"} (NR>1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%Tomato.bed}gat.out >${GAT}/${i%Tomato.bed}table.txt
awk -F"\t" 'BEGIN {OFS="\t"} (NR==1) {print FILENAME,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$9*-log($11)/log(10)}'  ${GAT}/${i%Tomato.bed}gat.out >${GAT}/${i%Tomato.bed}header.txt
cat ${GAT}/${i%Tomato.bed}header.txt ${GAT}/${i%Tomato.bed}table.txt > ${GAT}/${i%Tomato.bed}GAT_result_BL4+.txt
done
done

rm ${GAT}/*table.txt
rm ${GAT}/*header.txt
rm ${GAT}/*gat.out


cd -






