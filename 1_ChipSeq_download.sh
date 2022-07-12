###Author: Carina Gude
##download chromatin data - get the raw data from SRA
#Publication: preprint - "Dynamic epigenome changes in response to light in Brachypodium distachyon" ##NO: Fehlermeldung
#real Pubplication: https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.15288
#bsub -q bigmem -R "rusage[mem=20000]" -M 40000 /netscratch/dep_coupland/grp_turck/carina/Commands/Mtruncatula/1_Chip*.sh


cd /netscratch/dep_coupland/grp_turck/carina/Mtruncatula/

mkdir -p Downloads
cd Downloads 

#wget https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.8/Pecrix_et_al.ChIP-seq-peakcalling.zip
unzip Pecrix_et_al.ChIP-seq-peakcalling.zip
