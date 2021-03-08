# WEEK TWO

Zoom Meeting (2020.03.02 Tue): 
- Project aim: 
  - what are the main limitations in current scATAC-seq analysis?
  - for the reproducibility accessment of final peaks, what improvements does Chip-R have (or the adaptation of Chip-R on scATAC-seq data) on the analysis compared to other bioinformatic tools?
  - validation on robustness on clustering of cells
  - methods on cluster 1/2 (or 1/3, 1/4) of cells vs. all cells, analyze reproducibility scores also on the remaining cells
- Next week:
  - run scATAC-seq pipeline on raw data to process peaks
  - cell clustering and visualization
  - access reproducibility of peaks

### Reading List

- [Unsupervised clustering and epigenetic classification of single cells](https://www.nature.com/articles/s41467-018-04629-3)
- [Profiling chromatin regulatory landscape: insights into the development of ChIP-seq and ATAC-seq](https://link.springer.com/article/10.1186/s43556-020-00009-w)




### Limitations on current scATAC-seq processing
 - [Huidong.C et al 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5) scATAC-seq experiments sample DNA, which, due to low copy numbers (diploid in humans), lead to inherent data sparsity (1–10% of peaks detected per cell) compared to transcriptomic (scRNA-seq) data (10–45% of expressed genes detected per cell).
 - minor subpopulations are hard to detect unless deeper sequencing is performed.
 - complicated bioinformatic analysis on analyzing cell sequencing tracks individually and identify different clusters
 - going deeper into transcriptional mechanisms in individual cells will highlight potential therapeutic targets that are currently unknown
 - lack of open-source software for comprehensive processing, analysis, and visualization of scATAC data generated using all existing experimental protocols, users need to customize analysis step-by-step
 - there are only a limited set of tools for scATAC-seq data and a few reference datasets for cell type specific chromatin accessibility
 - [Feng.Y et al, 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3)  The chromatin accessibility at each base will be binary and the scATAC-seq data will be sparse because in diploid organisms, there are only two copies of DNA
 - [Mahdi.Z el al, 2018](https://www.nature.com/articles/s41467-018-04629-3) it is difficult to determine whether a region is absent in an individual cell due to the lack of openness or due to the sparse nature of data. This creates a challenging task in delineating distinct sub-populations, as only a few genomic regions will have overlapping reads in a large number of cells. If the population is unknown or marker genes are unavailable, then sub-population specific analysis becomes impractical
 - [Lei.X el al, 2019](https://www.nature.com/articles/s41467-019-12630-7) each open chromatin site of a diploid-genome single cell only has one or two opportunities to be captured. Normally, only a few thousand distinct reads (versus many thousands of possible open positions) are obtained per cell, thus resulting in many bona fide open chromatin sites of the cell that lack sequencing data signals (i.e., peaks). The analysis of scATAC-seq data hence suffers from the curse of “missingness” in addition to high dimensionality
   - chromVAR: only analyzes peaks in groups and lacks the resolution of individual peaks
   - scABC: heavily depends on landmark samples with high sequencing depths, and the Spearman rank can be ill-defined for data with many missing values
 -[Shaoqian.M et al, 2019](https://link.springer.com/article/10.1186/s43556-020-00009-w) not being able to identify rare peaks appearing only in scarce populations. One major limitation to current scATAC-seq approaches is that they capture only a subset of the open chromatin sites in single cells, a lot of sites may be lost or not detectable during both experimental and computational procedures




### Limitations of ChIP-R on scATAC-seq data





