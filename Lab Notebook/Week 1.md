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
- [snapshot](https://github.com/znavidi/scATAC-seq-analysis-pipeline)


**Analysis Pipeline:**

1. Data Processing:

  Download Data:
  
  ```
  wget https://cg.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_5k/atac_v1_pbmc_5k_fastqs.tar
  tar -cvf atac_v1_pbmc_5k_fastqs.tar  
  ```

  - output: atac_v1_pbmc_5k_S1_L001_I1_001.fastq.gz   
            atac_v1_pbmc_5k_S1_L001_R1_001.fastq.gz   
            atac_v1_pbmc_5k_S1_L001_R2_001.fastq.gz
            atac_v1_pbmc_5k_S1_L001_R3_001.fastq.gz   
            atac_v1_pbmc_5k_S1_L002_I1_001.fastq.gz   
            atac_v1_pbmc_5k_S1_L002_R1_001.fastq.gz
            atac_v1_pbmc_5k_S1_L002_R2_001.fastq.gz   
            atac_v1_pbmc_5k_S1_L002_R3_001.fastq.gz
  - I1 is the 8 bp sample barcode, 
  - R1 is the forward read, 
  - [R2 is the 16 bp 10x feature barcode](https://divingintogeneticsandgenomics.rbind.io/post/understand-10x-scrnaseq-and-scatac-fastqs/)
  - R3 is the reverse read
  - <ins>Question: Merge two lanes before QC/Reads Mapping or after?</ins>
  - <ins>Question: Need to add cell barcodes? (add R2 on R1 and R3)</ins>


  Attach cell barcodes:
  
  - input: fastq.gz
  - output: barcoded.fastq.gz
  - [sinto](https://timoast.github.io/sinto/index.html)

  
  ```
  sinto barcode -b 12 --barcode_fastq atac_v1_pbmc_5k_S1_L001_R2_001.fastq.gz \
  --read1 atac_v1_pbmc_5k_S1_L001_R1_001.fastq.gz \
  --read2 atac_v1_pbmc_5k_S1_L001_R3_001.fastq.gz 
  ```

  
  Quality Controll:

  - input: .fastq.gz
  - output: .trimmed.fastq.gz
  - FastQC, trimmomatic, and BWA-MEM

  ```
  ls *gz | while read id; do fastqc -o fastqc $id; done
  ```

  Reads Mapping:
    
  - input: .trammed.fastq.gz
  - output: .sam
  - reference genome: hg19
  - BWA MEM, Bowtie, STAR (faster, more accurate)
    
    
  ```
  # hg10 human reference genome
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  gzip -d hg19.fa.gz
  bwa index hg19.fa

  # reads maping

  ```
  
  Remove duplicated reads:
   
  - input:
  - output:
  - Picard

  ```
  
  ```

  Convert into BAM file, sort and index:
  
  - input: .sam
  - output: .bam
  - SAMTOOLS (sort, index, flagstat, view)
  - ?Generate and add unique barcodes for each single cell sample to sam files


  ```

  ```
  
   Create fragment file:
   
  - input: .bam
  - output: .fragments.bed
  - sinto fragments / cellranger-atac count

  ```

  ```
    
  Peak Calling:
    
  - input:
  - output:
  - Genrich / MACS2

  ```

  ```
    
  Access reproducibility of final peaks:
    
  - input: .bed
  - output: .bed
  - [Chip-R](https://github.com/rhysnewell/ChIP-R/)

  ```

  ```
   
  
2. Downstream Analysis:








**Dependencies**

Programming Languages:

- R
- Python

Software Packages:

- 


