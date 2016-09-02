# transcriptomics
Code used to conduct transcriptome analysis of Chalmydomonas reinhardtii and melanoma tumors
We used Trimmomatic to clean reads download from GEO (GSE34585). We used kallisto version 0.43.0 to obtain expression values. Custom Python scripts were used to call programs. 

To visualize the data, we used sklearn decomposition implementation in Python of PCA, and sklearn manifold implementation in Python of t-SNE. We used Rtsne implementation in R of t-SNE. 
