#!/bin/bash

cd ~

#variables and paths
path1=$(pwd)
mkdir -p ${path1}/GAT
OGAT=${path1}/GAT
GAT=${path1}/GAT
IDR=${path1}/IDR
TELOBED=${path1}/Telobox_BED
GENES=${path1}/Chrom_FeatureOverlap
ANNOTATION=/Annotation_Files
EXTRAS=/Extras
DOWNLOADS=/Downloads


#Find out annotations that are already in Genome GFF file#
awk '{print $3}'  ${DOWNLOADS}/*.gff | sort | uniq

cd ${ANNOTATION}
gff2bed < ${DOWNLOADS}/*gff > GFF3_annotation.bed

awk -F"\t" 'BEGIN {OFS="\t"} ($8=="gene") {print $0}' GFF3_annotation.bed > GFF3_genes.bed

##Annotation of TSS and TTS##
GENES=GFF3_genes.bed 
awk -vFS="\t" -vOFS="\t" '($6=="+")' ${GENES} > ${GENES}.for
awk -vFS="\t" -vOFS="\t" '($6=="-")' ${GENES} > ${GENES}.rev

TSS=${ANNOTATION}/TSS.bed
awk -vFS="\t" -vOFS="\t" '{ print $1, ($2-1), $2, $4, $5, $6, $7, $8, $9, $10 }' ${GENES}.for > ${TSS}.for
awk -vFS="\t" -vOFS="\t" '{ print $1, $3, ($3+1), $4, $5, $6, $7, $8, $9, $10 }' ${GENES}.rev > ${TSS}.rev
cat ${TSS}.for ${TSS}.rev | sort -k1,1 -k2,2n | sed 's/gene/TSS/g' > ${TSS}

TTS=${ANNOTATION}/TTS.bed
awk -vFS="\t" -vOFS="\t" '{ print $1, ($3-1), $3, $4, $5, $6, $7, $8, $9, $10 }' ${GENES}.for > ${TTS}.for
awk -vFS="\t" -vOFS="\t" '{ print $1, $2, ($2+1), $4, $5, $6, $7, $8, $9, $10 }' ${GENES}.rev > ${TTS}.rev
cat ${TTS}.for ${TTS}.rev | sort -k1,1 -k2,2n | sed 's/gene/TTS/g' > ${TTS}

##Annotation of promoters##
window=(3000 1000 100)
name=("3kbpromoter" "1kbpromoter" "100bpromoter")

for w in "${!window[@]}"
do
WINDOW="${window[w]}"
PROMOTERS="${name[w]}"
bedops --range -${WINDOW}:0 --everything ${TSS}.for > ${PROMOTERS}.for
bedops --range 0:${WINDOW} --everything ${TSS}.rev > ${PROMOTERS}.rev

#Reduce start location of promoters that overlap with a gene as e.g. the 3kb promoters could be jumping an entire gene
UPSTREAM_FILTERED_PROMOTERS="${name[w]}"_upstreamFiltered
bedmap --count --echo --echo-map-range ${PROMOTERS}.for ${GENES}.for | awk -vFS="|" -vOFS="\t" '{ if ($1==0) { print $2; } else { m=split($2,a,"\t"); tssStart=a[2]; tssEnd=a[3]; n=split($3,b,"\t"); geneEnd=b[3]; if ((geneEnd < tssEnd) && (geneEnd >= tssStart)) { print a[1], geneEnd, tssEnd, a[4], a[5], a[6], a[7], a[8], a[9], a[10]; } } }' > ${UPSTREAM_FILTERED_PROMOTERS}.for
bedmap --count --echo --echo-map-range ${PROMOTERS}.rev ${GENES}.rev | awk -vFS="|" -vOFS="\t" '{ if ($1==0) { print $2; } else { m=split($2,a,"\t"); tssStart=a[2]; tssEnd=a[3]; n=split($3,b,"\t"); geneStart=b[2]; if ((geneStart >= tssStart) && (geneStart < tssEnd)) { print a[1], tssStart, geneStart, a[4], a[5], a[6], a[7], a[8], a[9], a[10]; } } }' > ${UPSTREAM_FILTERED_PROMOTERS}.rev
bedops --everything ${UPSTREAM_FILTERED_PROMOTERS}.for ${UPSTREAM_FILTERED_PROMOTERS}.rev > ${UPSTREAM_FILTERED_PROMOTERS}.bed
done

##Annotation of Chromosomes##
#if chromosomes are annotated in GFF file:
awk -F"\t" 'BEGIN {OFS="\t"} ($8=="region") {print $0}' GFF3_annotation.bed | sed 's/region/chromosome/g' > GFF3_chromosomes.bed

#if chromosomes are not listed in GFF file: 
awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4="chromosome", $5=".", $6=".", $7="xxxx", $8="chromosome", $9=".", $10="ID=chromosome"}' ${EXTRAS}/CHROMLEN.bed > GFF3_chromosomes.bed


##Annotation of Exons and Introns##
awk -F"\t" 'BEGIN {OFS="\t"} ($8=="exon") {print $0}' GFF3_annotation.bed > GFF3_exons.bed
bedtools subtract -s -a GFF3_genes.bed -b GFF3_exons.bed | sed 's/gene/intron/g' > GFF3_introns.bed

for i in *promoter*upstreamFiltered.bed
do
bedtools subtract -A -s -a $i -b GFF3_genes.bed > GFF3_genes_${i%_upstreamFiltered.bed}.bed
done

#cut out the overlapping promoter pieces
bedtools subtract -s -a GFF3_genes_1kbpromoter.bed -b GFF3_genes_100bpromoter.bed | sed 's/gene/1kbpromoter/g' > GFF3_genes_1kbpromoter_noprox.bed
bedtools subtract -s -a GFF3_genes_3kbpromoter.bed -b GFF3_genes_1kbpromoter.bed | sed 's/gene/3kbpromoter/g' > GFF3_genes_3kbpromoter_no1kb.bed
sed -i 's/gene/100bpromoter/g' GFF3_genes_100bpromoter.bed

cat GFF3_genes_100bpromoter.bed GFF3_genes_1kbpromoter_noprox.bed GFF3_genes_3kbpromoter_no1kb.bed | sort -k1,1 -k2,2n >GFF3_promoter_unfiltered.bed

##Annotation of downstream region##
WINDOW=100
DOWNPROX=GFF3_genes_downprox.bed
bedops --range 0:${WINDOW} --everything ${TTS}.for > ${DOWNPROX}.for
bedops --range -${WINDOW}:0 --everything ${TTS}.rev > ${DOWNPROX}.rev
cat ${DOWNPROX}.for ${DOWNPROX}.rev | sort -k1,1 -k2,2n | sed 's/gene/downprox/g' > ${DOWNPROX}

##Annotation of Intergenic region##
FILTERED_DOWNPROX=GFF3_genes_downproxFiltered.bed
bedtools subtract -s -a ${DOWNPROX} -b GFF3_genes_100bpromoter.bed GFF3_genes.bed  > ${FILTERED_DOWNPROX}
bedtools subtract -s -a GFF3_promoter_unfiltered.bed -b ${FILTERED_DOWNPROX} > GFF3_promoter.bed
cat GFF3_promoter.bed GFF3_genes.bed GFF3_genes_downproxFiltered.bed | sort -k1,1 -k2,2n > GFF3_genes_Promoters_prox_down.bed
bedtools subtract -a ${ANNOTATION}/GFF3_chromosomes.bed -b GFF3_genes_Promoters_prox_down.bed | sed 's/chromosome/intergenic/g' > GFF3_intergenic.bed

#cat all files togehter and sort
cat GFF3_genes_Promoters_prox_down.bed GFF3_intergenic.bed | sort -k1,1 -k2,2n > GFF3_genes_promoters_intergenic.bed
cat GFF3_annotation.bed ${ANNOTATION}/TTS.bed ${ANNOTATION}/TSS.bed GFF3_genes_promoters_intergenic.bed GFF3_introns.bed | sort -k1,1 -k2,2n | awk -F '\t' '$2 != "-1"' > Mtruncatula_GAT.bed

rm -f *1kb* *100b* *.rev *.for  *3kb* *downprox* *+* *-* GFF3_[b-z]* TSS* TTS*
