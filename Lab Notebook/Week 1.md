# Week One

zoom meeting (2021.02.23 Tue) :

- Weekly meeting on Tuesday at 2pm Brisbane time
- project aims on "how peak-caller works on scATAC-seq data"
- general method on the project: define problems => find approaches to solve problems (make hypothesis) => do analysis


### Reading List:

- [ChIP-R: Assembling reproducible sets of ChIP-seq and ATAC-seq peaks from multiple replicates](https://www.biorxiv.org/content/10.1101/2020.11.24.396960v1.supplementary-material)

- [Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments](https://febs.onlinelibrary.wiley.com/doi/full/10.1016/j.febslet.2004.07.055)

- [Single-Cell Transcriptomic Analysis of Cardiac Differentiation from Human PSCs Reveals HOPX- Dependent Cardiomyocyte Maturation](https://linkinghub.elsevier.com/retrieve/pii/S1934590918304466)

- [Assessment of computational methods for the analysis of single-cell ATAC-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1854-5)

- [Profiling Chromatin Accessibility at Single-cell Resolution](https://www.sciencedirect.com/science/article/pii/S1672022921000115?via%3Dihub)

- [Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion. Nature Biotechnology](https://www.nature.com/articles/s41587-019-0206-z)



### **Dataset:**

[10X Genomics - Single cell ATAC](https://support.10xgenomics.com/single-cell-atac/datasets)

- [5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_pbmc_5k)    [summary](https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_5k/atac_v1_pbmc_5k_web_summary.html)
  - Peripheral blood mononuclear cells (PBMCs) from a healthy donor.
  - ~6700 transposed nuclei were loaded.
  - 4,654 nuclei were recovered.
  - Sequenced on Illumina NovaSeq with approximately 41k read pairs per cell.
  - 50bp read1, 8bp i7 (sample index), 16bp i5 (10x Barcode), 49bp read2.

### **scATAC-seq analysis tools:**

- [Signac](https://satijalab.org/signac/index.html)
- [MAESTRO](https://liulab-dfci.github.io/MAESTRO/example/ATAC_infrastructure_10x/ATAC_infrastructure_10x.html)
- [Sinto](https://timoast.github.io/sinto/index.html)
- [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
- [ArchR](https://www.archrproject.com/bookdown/getting-started-with-archr.html)
- [snapshot](https://github.com/znavidi/scATAC-seq-analysis-pipeline)
- [Scasat](https://github.com/ManchesterBioinference/Scasat)
- [chromVAR](https://greenleaflab.github.io/chromVAR/)
- [Brockman](https://carldeboer.github.io/brockman_pipe_example.html)
- [scATAC-pro](https://github.com/wbaopaul/scATAC-pro)
- [Cicero](https://cole-trapnell-lab.github.io/cicero-release/docs/)
- [Seurat](https://satijalab.org/seurat/index.html)

### **Analysis Pipeline:**

#### 1. Data Processing:

##### Download Data:
  
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

##### Attach cell barcodes:
  
  - input: fastq.gz
  - output: barcoded.fastq.gz
  - [sinto](https://timoast.github.io/sinto/index.html)
  > Storing the cell barcode in the read name is an easy way to track which reads came from which cells

  
  ```
  sinto barcode -b 12 --barcode_fastq atac_v1_pbmc_5k_S1_L001_R2_001.fastq.gz \
  --read1 atac_v1_pbmc_5k_S1_L001_R1_001.fastq.gz \
  --read2 atac_v1_pbmc_5k_S1_L001_R3_001.fastq.gz 
  
  cat atac_v1_pbmc_5k_S1_L001_R1_001.fastq.gz atac_v1_pbmc_5k_S1_L002_R1_001.fastq.gz > atac_v1_pbmc_5k_S1_R1_001.fastq.gz
  cat atac_v1_pbmc_5k_S1_L001_R3_001.fastq.gz atac_v1_pbmc_5k_S1_L002_R3_001.fastq.gz > atac_v1_pbmc_5k_S1_R3_001.fastq.gz
  ```

  
##### Quality Controll:

  - input: .fastq.gz
  - output: .trimmed.fastq.gz
  - FastQC, trimmomatic, and BWA-MEM

  ```
  ls *gz | while read id; do fastqc -o fastqc $id; done
 
  ```

##### Reads Mapping:
    
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
  bwa mem -t 8 ../hg19_reference/hg19.fa atac_v1_pbmc_5k_S1_R1_001.barcoded.fastq.gz atac_v1_pbmc_5k_S1_R3_001.barcoded.fastq.gz| samtools view -b - > atac_v1_pbmc_5k_S1_aln.bam
  ```
  
##### Remove duplicated reads:
   
  - input:
  - output:
  - Picard

  ```
  
  ```

##### Convert into BAM file, sort and index:
  
  - input: .sam
  - output: .bam
  - SAMTOOLS (sort, index, flagstat, view)
  - ?Generate and add unique barcodes for each single cell sample to sam files


  ```
  samtools sort -@ 12 atac_v1_pbmc_5k_S1_aln.bam -o atac_v1_pbmc_5k_S1_aln.sorted.bam
  samtools index -@ 12 atac_v1_pbmc_5k_S1_aln.sorted.bam
  ```
  
##### Create fragment file:
   
  - input: .bam
  - output: .fragments.bed
  - sinto fragments / cellranger-atac count

  ```
  sinto fragments -b atac_v1_pbmc_5k_S1_aln.sorted.bam -p 8 -f ../4_Fragment/atac_v1_pbmc_5k_S1_fragments.bed --barcode_regex "[^:]*"
  ```
    
##### Peak Calling:
    
  - input:
  - output:
  - Genrich / MACS2

  ```
  macs2 callpeak -t atac_v1_pbmc_5k_S1_L001_aln.sorted.bam -p 1e-3 --nomodel -n atac_v1_pbmc_5k_S1_L001_peak --outdir ../5_CallPeaks
  ```
   
  
#### 2. Downstream Analysis:

[Week3.md](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/Lab%20Notebook/Week%203.md)

