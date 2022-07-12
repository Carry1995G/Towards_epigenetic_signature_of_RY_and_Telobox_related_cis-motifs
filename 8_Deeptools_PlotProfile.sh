#!/bin/bash
##Generate Plot Profiles from Matrix 
#bsub -q normal -R "rusage[mem=10240]" -M 12288 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/8_Deeptools_Plot*

cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

mkdir -p Plots
cd Matrix
echo "Generate plot profile" 

##quickly generate plot profiles without labels 
#for i in *.gz
#do
#plotProfile -m $i  \
#			-o ../Plots/${i%.gz}_Profileplot.png \
#			--samplesLabel "${i%_Tomato.gz} motifs in Solanum lycopersicum"
#done

##quickly generate plot profiles with labels
#plotProfile -m *.gz \
#			-o ../Plots/Telobox+Telolike_Profileplot.png \
#			--plotTitle "Motifs in Solanum lycopersicum" \
#			--samplesLabel "Telobox motifs" "Telolike motifs"

##With Blacklists###

for i in *_Mtruncatula.gz
do
plotProfile -m $i  \
			-o ../Plots/${i%_Mtruncatula.gz}_Profileplot_$p.png \
			--samplesLabel "${i%_Mtruncatula.gz} motifs" \
			--plotType se \
			-y "Coverage per 50 bp" \
			--legendLocation "best"
done


#for p in ${plotType[@]}
#do
#plotProfile -m Motifs_Mtruncatula.gz \
#			-o ../Plots/Motifs_Profileplot_$p.png \
#			--plotTitle "Motifs in Medicago truncatula " \
#			--samplesLabel "RY motifs" "Telobox motifs" "Telolike motifs" \
#			--plotType se
#done

#for i in *Telobox_Telolike*BL.gz
#do
#plotProfile -m $i \
#			-o ../Plots/${i%.gz}_Profileplot.png \
#			--plotTitle "Motifs in Solanum lycopersicum" \
#			--samplesLabel "Telobox motifs" "Telolike motifs"
#done