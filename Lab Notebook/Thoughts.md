**COURSE**

- [ ] I've received the Electronic Course Profile and noted the important due dates. 
- [ ] Will all the assessment submitted through the blackboard? (The blackboard has not updated for this research course yet)
- [ ] the project plan is due on about 5th week, so should I get the data and start programming before or after due date of project plan？
- [ ] there is a 20% "Laboratory Notebook and assessment of research performance", what is the best way to keep my daily record so the examiners and you can see?
- [ ] So the presentation and report was due on June 7th, does it mean this research will end on June 7th?

**RSEARCH**

- [ ] when will I get the scATAC dataset? (Is it human or other species?)
- [ ] what programming language will be used most often?
- [ ] Do I work on the same/different hpc server? s4606549@bioinfr.server.science.uq.edu.au



- [ ] 



***

[TOC]

[Github](https://github.com/Fiona-Pan/Master-Research-Project)

***

## **ATAC-seq**

**chromatin accessibility is dynamic**

- Closed chromatin: inactive genes and enhancers
- open chromatin: active genes and enhancers
  - These open regions can be sequenced selectively
  - open chromatin can:
    - quickly determine if epigenetics is involved in response
    - Infer which histone modification(s) are of interest based on correlation of open chromatin regions to genomic annotations
    - Identify what transcription factors are driving cell fate, disease or response
    - Machanistic insight into how your pathway/cellular response/disease is regulated

##### LIMITATIONS

- Standard ATAC-seq averages out cellular heterogeneity such that signals from rare populations or important cell subtypes are masked
- challenge to collect rare cell populations to achieve the recommended numbers for ATAC-seq

##### STEPS

1. Tn5 transposase gains access to open/accessible regions of the genome
2. Transposase inserts known DNA sequences (adapters) into open regions
3. Libraries are amplified and sequenced<img src="/Users/apple/Library/Application Support/typora-user-images/image-20210218115350874.png" alt="image-20210218115350874" style="zoom:25%;" /><img src="/Users/apple/Library/Application Support/typora-user-images/image-20210218115440624.png" alt="image-20210218115440624" style="zoom: 25%;" />

##### SUMMARY

- ATAC-seq map open chromatin on a genome-wide level
- ATAC-seq allows for the interrogation of open and accessible regions of the genome
- ATAC-seq reveals all active loci (promoters and enhancers) genome wide.
- ATAC-seq data provide distinct peaks representing specific active loci
- The ease of the assay allows for rapid detection of aberrantly active loci
- ATAC-seq can detect active transcription factor binding sites within promoters, informing choice of drug candidates

***

## **Single-Cell ATAC-Seq**

- scATAC-seq can deconvolute heterogeneity of a mixed population
- utilizes Tn5 transposase and barcoding of individual cells to profile chromatin accessibility at single cell resolution
- single cells can be captured by combinatorial cell indexing strategies or with the use of a microfluidic device
- scATAC-seq data is represented using tSNE plots<img src="/Users/apple/Library/Application Support/typora-user-images/image-20210218121318758.png" alt="image-20210218121318758" style="zoom:25%;" />
- scATAC-seq can help address many experimental questions
  - classify different cell populations in a heterogenous population
  - understand the tumor microenvironment
  - examine cell differentiation and cell trajectories during development
  - Identify gene regulatory networks involved in development or disease
- <img src="/Users/apple/Library/Application Support/typora-user-images/image-20210218121710167.png" alt="image-20210218121710167" style="zoom:50%;" />

***

# PAPERS

#### [ChIP-R: Assembling reproducible sets of ChIP-seq and ATAC-seq peaks from multiple replicates](https://www.biorxiv.org/content/10.1101/2020.11.24.396960v1.supplementary-material)

[GITHUB LINK](https://github.com/rhysnewell/ChIP-R/)

- chip-seq detect genome-wide DNA-protein interactions, a key tool for understanding transcriptional regulation
  - limitations: 
    - low specificity of antibody & cellular heterogeneity of sample => cause “peak” callers to output noise and experimental artefacts
    - ChIP-seq protocols are susceptible to produce false positive signals as DNA regions not bound by the target protein can be pulled down indiscriminately during immunoprecipitation
- Statistically combining multiple experimental replicates from the same condition could signifi- cantly enhance our ability to distinguish actual transcription factor binding events, even when peak caller accuracy and consistency of detection are compromised.
- rank-product test => statistically evaluate the reproducibility from any number of ChIP-seq experimental replicates
- "ChIP-R" extends to evaluate ATAC-seq peaks => finding reproducible peak sets even at low sequencing depth





## **ATAC-seq**

#### [Single-Cell Transcriptomic Analysis of Cardiac Differentiation from Human PSCs Reveals HOPX- Dependent Cardiomyocyte Maturation](https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(18)30446-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1934590918304466%3Fshowall%3Dtrue)

- performed exten- sive single-cell transcriptomic analyses to map fate choices and gene expression programs during car- diac differentiation of hPSCs and identified strate- gies to improve in vitro cardiomyocyte differentiation.
- identify the non- DNA binding homeodomain protein *HOPX*, a key regulator of heart development and hypertrophy as dysregu- lated during differentiation and a potential cause for the imma- ture state of hPSC-derived cardiomyocytes *in vitro*





## **RANK-PRODUCT**

[RankProd R Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/RankProd.html) 		[Package Info PDF](http://127.0.0.1:15607/library/RankProd/doc/RankProd.pdf)

#### [Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments](https://febs.onlinelibrary.wiley.com/doi/full/10.1016/j.febslet.2004.07.055)

[link](http://home.cc.umanitoba.ca/~psgendb/birchhomedir/doc/MeV/manual/rp.html)

- identifying differentially expressed genes that does not originate from a sophisticated statistical model but rather from an analysis of biological reasoning
- based on calculating rank products (RP) from replicate experiments
- 

















