# Towards epigenetic signature of RY and Telobox related cis-motifs 

Welcome,

here you can find the code I more detailed description on the code I used for my master's thesis.


## Workflow:

1. [Data Download & Blacklists](/1.Data_Download_Blacklists.bash)

These two scripts can be run in parallel:

2. [Cis-motif Annotation](/2a.Cis-motif_data_preparation.bash)  OR
2. [Chromatin data preparation](/2b.Chromatin_data_preparation.sh)

As soon as 2a and 2b are done, this script can be used:
3. Motif distribution in H3K27me3+ and H3K27me3- genes 

Next comes the GAT Analysis with 4.3 run in R.
4. GAT-Annotation
5. GAT-Analysis 
6. GAT-Visualization

In the end comes the GO-Analysis with 5.2 and 5.3 run in R.
7. GO-Preparation
8. TopGO-Analysis
9. GO visualization
