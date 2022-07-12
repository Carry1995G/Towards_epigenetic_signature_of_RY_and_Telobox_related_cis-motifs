#!/bin/bash
#bsub -q normal -R "rusage[mem=10240]" -M 12288 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/11_GAT_Annotation.sh

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

####################################################################
#variables and paths
path1=$(pwd)
mkdir -p ${path1}/GAT
OGAT=${path1}/GAT
GAT=${path1}/GAT
IDR=${path1}/IDR
TELOBED=${path1}/Telobox_BED
GENES=${path1}/Chrom_FeatureOverlap
ANNOTATION=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Annotation_Files
EXTRAS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Extras
DOWNLOADS=/netscratch/dep_coupland/grp_turck/carina/Mtruncatula/Downloads


####Annotation already in file####
#awk '{print $3}' *1.8.gff3 | sort | uniq
#CDS
#LTR_retrotransposon
#RR_tract
#TRIM_retrotransposon
#exon
#five_prime_UTR
#gene
#genes
#inverted_repeat
#long_terminal_repeat
#mRNA
#miRNA
#ncRNA
#pre_miRNA
#primer_binding_site
#protein_match
#rRNA
#repeat_region
#tRNA
#target_site_duplication
#terminal_inverted_repeat
#terminal_inverted_repeat_element
#three_prime_UTR

#Annotation of promoters
cd ${ANNOTATION}
gff2bed < ${DOWNLOADS}/*1.8.gff3 > GFF3_annotation.bed

awk -F"\t" 'BEGIN {OFS="\t"} ($8=="gene") {print $0}' GFF3_annotation.bed > GFF3_genes.bed

GENES=GFF3_genes.bed #Mtruncatula.bed
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

window=(3000 1000 100)
name=("3kbpromoter" "1kbpromoter" "100bpromoter")

for w in "${!window[@]}"
do
WINDOW="${window[w]}"
PROMOTERS="${name[w]}"
bedops --range -${WINDOW}:0 --everything ${TSS}.for > ${PROMOTERS}.for
bedops --range 0:${WINDOW} --everything ${TSS}.rev > ${PROMOTERS}.rev


UPSTREAM_FILTERED_PROMOTERS="${name[w]}"_upstreamFiltered
bedmap --count --echo --echo-map-range ${PROMOTERS}.for ${GENES}.for | awk -vFS="|" -vOFS="\t" '{ if ($1==0) { print $2; } else { m=split($2,a,"\t"); tssStart=a[2]; tssEnd=a[3]; n=split($3,b,"\t"); geneEnd=b[3]; if ((geneEnd < tssEnd) && (geneEnd >= tssStart)) { print a[1], geneEnd, tssEnd, a[4], a[5], a[6], a[7], a[8], a[9], a[10]; } } }' > ${UPSTREAM_FILTERED_PROMOTERS}.for
bedmap --count --echo --echo-map-range ${PROMOTERS}.rev ${GENES}.rev | awk -vFS="|" -vOFS="\t" '{ if ($1==0) { print $2; } else { m=split($2,a,"\t"); tssStart=a[2]; tssEnd=a[3]; n=split($3,b,"\t"); geneStart=b[2]; if ((geneStart >= tssStart) && (geneStart < tssEnd)) { print a[1], tssStart, geneStart, a[4], a[5], a[6], a[7], a[8], a[9], a[10]; } } }' > ${UPSTREAM_FILTERED_PROMOTERS}.rev

bedops --everything ${UPSTREAM_FILTERED_PROMOTERS}.for ${UPSTREAM_FILTERED_PROMOTERS}.rev > ${UPSTREAM_FILTERED_PROMOTERS}.bed
done

######### file 2 #######

# I used a scipt on the internet to define the promotors but not all of them are fully in intergenic space, in particular the 3kb promoters may be jumping an entire gene!
cd ${ANNOTATION}
awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4="chromosome", $5=".", $6=".", $7="maker_ITAG", $8="chromosome", $9=".", $10="ID=chromosome"}' ${EXTRAS}/Mtruncatula_chromlen.bed > GFF3_chromosomes.bed

#awk -F"\t" 'BEGIN {OFS="\t"} ($8=="region") {print $0}' GFF3_annotation.bed | sed 's/region/chromosome/g' > GFF3_chromosomes.bed
awk -F"\t" 'BEGIN {OFS="\t"} ($8=="exon") {print $0}' GFF3_annotation.bed > GFF3_exons.bed
bedtools subtract -s -a GFF3_genes.bed -b GFF3_exons.bed | sed 's/gene/intron/g' > GFF3_introns.bed


for i in *promoter*upstreamFiltered.bed
do
bedtools subtract -A -s -a $i -b GFF3_genes.bed > GFF3_genes_${i%_upstreamFiltered.bed}.bed
done

#now it makes sense to cut out the overlapping pieces

bedtools subtract -s -a GFF3_genes_1kbpromoter.bed -b GFF3_genes_100bpromoter.bed | sed 's/gene/1kbpromoter/g' > GFF3_genes_1kbpromoter_noprox.bed
bedtools subtract -s -a GFF3_genes_3kbpromoter.bed -b GFF3_genes_1kbpromoter.bed | sed 's/gene/3kbpromoter/g' > GFF3_genes_3kbpromoter_no1kb.bed
sed -i 's/gene/100bpromoter/g' GFF3_genes_100bpromoter.bed

cat GFF3_genes_100bpromoter.bed GFF3_genes_1kbpromoter_noprox.bed GFF3_genes_3kbpromoter_no1kb.bed | sort -k1,1 -k2,2n >GFF3_promoter_unfiltered.bed


WINDOW=100
DOWNPROX=GFF3_genes_downprox.bed
bedops --range 0:${WINDOW} --everything ${TTS}.for > ${DOWNPROX}.for
bedops --range -${WINDOW}:0 --everything ${TTS}.rev > ${DOWNPROX}.rev
cat ${DOWNPROX}.for ${DOWNPROX}.rev | sort -k1,1 -k2,2n | sed 's/gene/downprox/g' > ${DOWNPROX}

FILTERED_DOWNPROX=GFF3_genes_downproxFiltered.bed
bedtools subtract -s -a ${DOWNPROX} -b GFF3_genes_100bpromoter.bed GFF3_genes.bed  > ${FILTERED_DOWNPROX}
bedtools subtract -s -a GFF3_promoter_unfiltered.bed -b ${FILTERED_DOWNPROX} > GFF3_promoter.bed

cat GFF3_promoter.bed GFF3_genes.bed GFF3_genes_downproxFiltered.bed | sort -k1,1 -k2,2n > GFF3_genes_Promoters_prox_down.bed


#bedtools subtract -s -a GFF3_chromosomes.bed -b GFF3_genes_Promoters_prox_down.bed | sed 's/chromosome/intergenic/g' > GFF3_intergenic.bed
#awk -F"\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6="+",$7,$8,$9,$10}' GFF3_chromosomes.bed > chromosomes+.bed
#awk -F"\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6="-",$7,$8,$9,$10}' GFF3_chromosomes.bed > chromosomes-.bed
bedtools subtract -a ${ANNOTATION}/GFF3_chromosomes.bed -b GFF3_genes_Promoters_prox_down.bed | sed 's/chromosome/intergenic/g' > GFF3_intergenic.bed
#bedtools subtract -s -a ${ANNOTATION}/chromosomes-.bed -b GFF3_genes_Promoters_prox_down.bed | sed 's/chromosome/intergenic/g' > GFF3_-intergenic.bed
#cat GFF3_+intergenic.bed GFF3_-intergenic.bed | sort -k1,1 -k2,2n > GFF3_intergenic.bed

cat GFF3_genes_Promoters_prox_down.bed GFF3_intergenic.bed | sort -k1,1 -k2,2n > GFF3_genes_promoters_intergenic.bed
cat GFF3_annotation.bed ${ANNOTATION}/TTS.bed ${ANNOTATION}/TSS.bed GFF3_genes_promoters_intergenic.bed GFF3_introns.bed | sort -k1,1 -k2,2n | awk -F '\t' '$2 != "-1"' > Mtruncatula_GAT.bed

rm -f *1kb* *100b* *.rev *.for  *3kb* *downprox* *+* *-* GFF3_[b-z]* TSS* TTS*

bsub -q normal -R "rusage[mem=10240]" -M 12288 bash /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/12_GAT.sh