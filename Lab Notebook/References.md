# Reading List

[Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments](https://febs.onlinelibrary.wiley.com/doi/full/10.1016/j.febslet.2004.07.055)

- identifying differentially expressed genes that does not originate from a sophisticated statistical model but rather from an analysis of biological reasoning
- based on calculating rank products (RP) from replicate experiments


[Assessment of computational methods for the analysis of single-cell ATAC-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5)
- [Github](https://github.com/pinellolab/scATAC-benchmarking/)

[Profiling Chromatin Accessibility at Single-cell Resolution](https://www.sciencedirect.com/science/article/pii/S1672022921000115?via%3Dihub)
- scATAC-seq processing constitutes the following five sequential steps
1. Processing raw sequencing data
 - trimming sequencing adapters, 
 - eliminating poor quality reads, 
 - mapping paired-end reads, 
 - eliminating cells with a library size that falls below a chosen threshold, 
 - aggregating cell barcodes for downstream processing
2. Pre-defined feature selection
 - the overall variance in a dataset is reduced by restricting the analysis to pre-defined motifs or a list of genomic/gene set annotations
   - [chromVAR](http://www.github.com/GreenleafLab/chromVAR) (tool): aggregated reads based on known TF motifs
   - converted sparse accessibility by cell barcode matrix to a more stable, bias-corrected deviation per TF motif matrix, enabling cell-to-cell similarity measurements by analyzing gain or loss of accessibility within defined genomic features
3. Heterogeneity calculation of single cells
 - Using the filtered accessibility–barcode matrix, a dissimilarity measurement that quantitates the extent of divergence between two cells is calculated based on differential peak accessibility.
 - [SCRAT](https://github.com/zji90/SCRAT) (tool): first aggregates co-activated sites to derive pathway level accessibility. It then uses pathway-level accessibility as a feature to compute cell dissimilarity
 - It is shown that aggregation can increase signal-to-noise ratio and lead to improved dissimilarity measures that better separate cells by cell type
4. Dimensionality reduction
 - reducing its dimensionality, such that biological variance is retained but random variables contributing to noise are eliminated, is critical for downstream analysis and visualization
   - MDS
   - tSNE
   - UMAP
5. Cell clustering and differential accessibility analysis
 -  hierarchically cluster cells into different groups
 -  A number of different statistical tests (i.e., Fisher exact, Information gain) can then be employed to calculate differential accessibility signatures between two clusters.
 -  (1) mapping co-accessible DNA elements to connect regulatory sequences with its targets, 
 -  (2) predicting cell-type specific transcription factor activity, and 
 -  (3) integrating scATAC-seq data with scRNA-seq, spatial transcriptomics, and other single-cell measurements





