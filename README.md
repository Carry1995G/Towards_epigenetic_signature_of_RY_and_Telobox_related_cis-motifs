# Towards epigenetic signature of RY and Telobox related cis-motifs 

Welcome,

here you can find the code I more detailed description on the code I used for my master's thesis.


## Workflow:

1. [Data Download & Blacklists](/1.Data_Download_Blacklists.bash)


These scripts can be run in parallel:

2a Cis-motif Annotation   OR

2b Chromatin data preparation

As soon as 2a and 2b are done, this script can be used:

3. Motif distribution in H3K27me3+ and H3K27me3- genes 

Next comes the GAT Analysis with 4.3 run in R.
4.1 GAT-Annotation
4.2 GAT-Analysis 
4.3 GAT-Visualization

In the end comes the GO-Analysis with 5.2 and 5.3 run in R.
5.1 GO-Preparation
5.2 TopGO-Analysis
5.3 GO visualization
