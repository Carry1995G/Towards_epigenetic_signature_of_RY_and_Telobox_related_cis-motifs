# Towards epigenetic signature of RY and Telobox related cis-motifs 

Welcome,

here you can find the more detailed description on the code I used for my master's thesis.


## Workflow:

1. [Data Download & Blacklists](/1.Data_Download_Blacklists.bash)

These two scripts can be run in parallel:

2. [Cis-motif Annotation](/2a.Cis-motif_data_preparation.bash)  OR
2. [Chromatin data preparation](/2b.Chromatin_data_preparation.sh)

This script can be used after the first three scripts have run.

3. [Motif distribution in H3K27me3+ and H3K27me3- genes](/3.Motif_distribution_in_H3K27me3+_and_H3K27me3-_genes.sh)

Next comes the GAT Analysis with 6 run in R.

4. [GAT-Annotation](/11_GAT_Annotation.sh)
5. [GAT-Analysis](/12_GAT.sh)
6. [GAT-Visualization](/7.GAT_Visualization.R)

In the end comes the GO-Analysis with 8 and 9 run in R.

7. [GO-Preparation](/8_GO_Preparation.sh)
8. [TopGO-Analysis](/9.TopGO-analysis.R)
9. [GO visualization](/GO-analysis_visualization.R)
