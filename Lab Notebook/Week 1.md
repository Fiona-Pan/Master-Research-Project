# Week One


Reading List:

- [ChIP-R: Assembling reproducible sets of ChIP-seq and ATAC-seq peaks from multiple replicates](https://www.biorxiv.org/content/10.1101/2020.11.24.396960v1.supplementary-material)

- [Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments](https://febs.onlinelibrary.wiley.com/doi/full/10.1016/j.febslet.2004.07.055)

- [Single-Cell Transcriptomic Analysis of Cardiac Differentiation from Human PSCs Reveals HOPX- Dependent Cardiomyocyte Maturation](https://linkinghub.elsevier.com/retrieve/pii/S1934590918304466)

- [Assessment of computational methods for the analysis of single-cell ATAC-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5)

- [Profiling Chromatin Accessibility at Single-cell Resolution](https://www.sciencedirect.com/science/article/pii/S1672022921000115?via%3Dihub)

- [Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion. Nature Biotechnology](https://www.nature.com/articles/s41587-019-0206-z)

zoom meeting with Mikael (2021.02.23 Tue) :

- Weekly meeting on Tuesday at 2pm Brisbane time
- project aims on "how peak-caller works on scATAC-seq data"
- general method on the project: define problems => find approaches to solve problems (make hypothesis) => do analysis


**Dataset:**

- public dataset from 10X Genomics: [5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_pbmc_5k)

**scATAC-seq analysis tools:**

- [Signac](https://satijalab.org/signac/index.html)
- [MAESTRO](https://liulab-dfci.github.io/MAESTRO/example/ATAC_infrastructure_10x/ATAC_infrastructure_10x.html)
- [Sinto](https://timoast.github.io/sinto/index.html)
- [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
- [ArchR](https://www.archrproject.com/bookdown/getting-started-with-archr.html)


**Analysis Pipeline:**

1. Data Processing:

- FASTQs &rarr; Trimming & Filtering &rarr; BAM, fragments.txt

  Download Data:
  ```
  wget https://cg.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_5k/atac_v1_pbmc_5k_fastqs.tar
  tar -cvf atac_v1_pbmc_5k_fastqs.tar  
  ```
  --create (-c) - Create a new tar archive. 
  
  --verbose (-v) - Show the files being processed by the tar command. 
  
  --file=archive=name (-f archive-name) - Specifies the archive file name.

2. Downstream Analysis:


**Dependencies**

Programming Languages:

- R
- Python

Software Packages:

- 


