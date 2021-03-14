# Reading List

- [Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments](https://febs.onlinelibrary.wiley.com/doi/full/10.1016/j.febslet.2004.07.055)
  - calculating rank products (RP) from replicate experiments
  - identifying differentially expressed genes that does not originate from a sophisticated statistical model but rather from an analysis of biological reasoning
  - provides a straightforward and statistically stringent way to determine the significance level for each gene and allows for the flexible control of the false-detection rate and familywise error rate in the multiple testing situation of a microarray experiment

- [Assessment of computational methods for the analysis of single-cell ATAC-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5)
  - [Github](https://github.com/pinellolab/scATAC-benchmarking/)
  - computational challenges associated with scATAC-seq analysis: inherently sparse data, determination of features, peak calling, the effects of sequencing coverage and noise, clustering performance
  - analysis steps: 
   - raw reads obtained in .fastq format for each single cell are mapped to a reference genome, producing aligned reads in .bam format
   - peak calling and read counting return the genomic position and the read counts in .bed and .txt format
   - data in these file formats in then used for downstream analysis
   - feature matrix constructed from single cells (by counting the number of reads at each peak for every cell) followed by visualization, clustering, trrajectory inference, determination of differential accessibility, prediction of cis-regulatory networks
  - public data with aligned files in BAM format, using three commonly used clustering approaches (K-means, Louvain, hierarchical clustering) and UMAP projection to find putative subpopulations and visualize cell-to-cell similarities for each method
  - computational methods include: BROCKMAN, chromVAR, Cicero, cisTopic, *Cusanovich2018*, Gene Scoring, scABC, Scasat, SCRAT, SnapATAC
   - single-BAM file &rarr; feature matrix construction (define regions, count features, transformation, dimensionality reduction)

- [Profiling Chromatin Accessibility at Single-cell Resolution](https://www.sciencedirect.com/science/article/pii/S1672022921000115?via%3Dihub)

- [Mapping genetic effects on cell type-specific chromatin accessibility and annotating complex trait variants using single nucleus ATAC-seq](https://www.biorxiv.org/content/10.1101/2020.12.03.387894v1)

- [Single-Cell Transcriptomic Analysis of Cardiac Differentiation from Human PSCs Reveals HOPX- Dependent Cardiomyocyte Maturation](https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(18)30446-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1934590918304466%3Fshowall%3Dtrue)

- [ChIP-R: Assembling reproducible sets of ChIP-seq and ATAC-seq peaks from multiple replicates](https://www.biorxiv.org/content/10.1101/2020.11.24.396960v1.supplementary-material)
  - [Github](https://github.com/rhysnewell/ChIP-R/)

- [Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion. Nature Biotechnology](https://www.nature.com/articles/s41587-019-0206-z)

- [Unsupervised clustering and epigenetic classification of single cells](https://www.nature.com/articles/s41467-018-04629-3)

- [Profiling chromatin regulatory landscape: insights into the development of ChIP-seq and ATAC-seq](https://link.springer.com/article/10.1186/s43556-020-00009-w)

- [Single-cell ATAC sequencing analysis: From data preprocessing to hypothesis generation](https://www.sciencedirect.com/science/article/pii/S2001037020303019?via%3Dihub)

- [From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3)

- [Single-cell ATAC sequencing analysis: From data preprocessing to hypothesis generation](https://www.sciencedirect.com/science/article/pii/S2001037020303019)





