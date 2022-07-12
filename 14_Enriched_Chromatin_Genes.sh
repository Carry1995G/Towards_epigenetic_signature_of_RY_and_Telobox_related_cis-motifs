#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/14_Enriched_Chromatin*

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

mkdir -p Chrom_FeatureOverlap
path1=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula
Downloads=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Downloads
Genes=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Chrom_FeatureOverlap
Annot=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Annotation_Files
#Blacklists=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Blacklists

mkdir -p Chrom_FeatureOverlap/Genes
GENE=${Genes}/Genes

cd ${Downloads}
########### GENES#########
#Mtruncatula.bed i
gff2bed < ${Downloads}/Mtrun*1.8.gff3 | awk '$8=="gene" {print $0}'> ${Annot}/GFF2BED_genes.bed


#Annotationfile im BED format# Hier nicht nötig, da genes im Text vorkommen 
#awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";";}' *gtf |\
#	convert2bed --input=gtf - > MtruncatulaBEDAnnotation.bed

cd ${path1}/Peaks

for i in *intersect.bed
do
bedtools intersect -a ${Annot}/GFF2BED_genes.bed -b $i -f 0.5 -u > ${Genes}/${i%_SICER_intersect.bed}K27+Genes.bed
awk -vFS="\t" -vOFS="\t" 'split($4, a, /[:]/) >= 2 {print a[2]}' ${Genes}/${i%_SICER_intersect.bed}K27+Genes.bed | \
sed 's/MtrunA17/&_/' - | awk /MtrunA17_Chr[1-9]/ | uniq > ${GENE}/${i%_SICER_intersect.bed}K27+Genes.txt
done

for i in *intersect.bed
do
bedtools intersect -a ${Annot}/GFF2BED_genes.bed -b $i -f 0.5 -v > ${Genes}/${i%_SICER_intersect.bed}K27-Genes.bed
grep -w 'gene' ${Genes}/${i%_SICER_intersect.bed}K27-Genes.bed > ${Genes}/${i%_SICER_intersect.bed}.temp
perl -ane 'print /ID=gene:"([A-Za-z0-9]+)"Alias=/, "\n";' ${Genes}/${i%_SICER_intersect.bed}.temp | uniq > ${GENE}/${i%_SICER_intersect.bed}K27-Genes.txt
done

#########CDS########
cd ${Downloads}
gff2bed < ${Downloads}/Mtrun*1.8.gff3 | awk '$8=="CDS" {print $0}'> ${Annot}/GFF2BED_CDS.bed

cd ${path1}/Peaks

bedtools intersect -a ${Annot}/GFF2BED_CDS.bed -b *intersect.bed -f 0.5 -u > ${Genes}/H3K27me3+CDS.bed
grep -w 'CDS' ${Genes}/H3K27me3+CDS.bed > ${Genes}/H3K27me3+CDS.temp
perl -ane 'print /ID=gene:"([A-Za-z0-9]+)";/, "\n";' ${Genes}/H3K27me3+CDS.temp | uniq > ${GENE}/H3K27me3+CDS.txt


bedtools intersect -a ${Annot}/GFF2BED_CDS.bed -b *intersect.bed -f 0.5 -v > ${Genes}/H3K27me3-CDS.bed
grep -w 'CDS' ${Genes}/H3K27me3-CDS.bed > ${Genes}/H3K27me3-CDS.temp
perl -ane 'print /ID=gene:"([A-Za-z0-9]+)"Alias=/, "\n";' ${Genes}/H3K27me3-CDS.temp | uniq > ${GENE}/H3K27me3-CDS.txt



rm ${Genes}/*temp
bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/15_Deeptools_He*


#########################################################################################
#for i in Annotate*.txt
#do
#cut -f4 $i | grep -o '\(A[TG0-9]*\)' | cut -c-9 | sed "1 d" |sort | uniq >$i.AGIs.txt
#done

#for i in Annotate*FDR.01.txt
#do 
#awk -F"\t" 'BEGIN {OFS="\t"}{if (/tss/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"TSS",$4;
#	else if (/5*UTR/)  print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"5UTR",$4;
#	else if (/3*UTR/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"3UTR",$4 ;
#	else if (/TTS/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"TTS",$4 ;
#	else if (/exon/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"genebody",$4;
#	else if (/intron/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"genebody",$4 ;
#	else if (/1kb-promoter/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"1kb-promoter",$4 ;
#	else if (/3kb-promoter/) print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"3kb-promoter",$4 ;
#	else print "chr"$1,$2,$3,substr( $4, index($4,"A"), 9 ),"Intergenic",$4 ;}' $i >$i.final


#wc -l $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="genebody" { count++ } END { print "genebody" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="3UTR" { count++ } END { print "3UTR" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="TSS" { count++ } END { print "TSS" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="5UTR" { count++ } END { print "5UTR" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="TTS" { count++ } END { print "TTS" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="1kb-promoter" { count++ } END { print "1kb-promoter" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="3kb-promoter" { count++ } END { print "3kb-promoter" , count }' $i.final >>Homer.analysis.txt
#awk -F"\t" '$5=="Intergenic" { count++ } END { print "Intergenic" , count }' $i.final >>Homer.analysis.txt
#done