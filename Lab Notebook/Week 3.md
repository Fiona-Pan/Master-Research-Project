# WEEK 3


Zoom meeting (2020.03.09 Tue):

- read pre-processed dataset in R (peak-by-cell matrix)
- clustering and find biological meanings

### Snaptools Pipeline

- input: raw fastq files
- output: peak-by-cell matrix
- [Github](https://github.com/r3fang/SnapATAC)

```
# download data
wget https://cg.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_5k/atac_v1_pbmc_5k_fastqs.tar
tar -cvf atac_v1_pbmc_5k_fastqs.tar 

# barcode demultipleing
snaptools dex-fastq \
--input-fastq=atac_v1_pbmc_5k_S1_L001_R1_001.fastq.gz \
--output-fastq=atac_v1_pbmc_5k_S1_L001_R1_001.dex.fastq.gz \
--index-fastq-list atac_v1_pbmc_5k_S1_L001_R2_001.fastq.gz

snaptools dex-fastq --input-fastq=atac_v1_pbmc_5k_S1_L002_R1_001.fastq.gz \
--output-fastq=atac_v1_pbmc_5k_S1_L002_R1_001.dex.fastq.gz \
--index-fastq-list atac_v1_pbmc_5k_S1_L002_R2_001.fastq.gz

snaptools dex-fastq --input-fastq=atac_v1_pbmc_5k_S1_L001_R3_001.fastq.gz \
--output-fastq=atac_v1_pbmc_5k_S1_L001_R3_001.dex.fastq.gz \
--index-fastq-list atac_v1_pbmc_5k_S1_L001_R2_001.fastq.gz

snaptools dex-fastq --input-fastq=atac_v1_pbmc_5k_S1_L002_R3_001.fastq.gz \
--output-fastq=atac_v1_pbmc_5k_S1_L002_R3_001.dex.fastq.gz \
--index-fastq-list atac_v1_pbmc_5k_S1_L002_R2_001.fastq.gz

# combine lanes
cat atac_v1_pbmc_5k_S1_L001_R1_001.dex.fastq.gz atac_v1_pbmc_5k_S1_L002_R1_001.dex.fastq.gz > atac_v1_pbmc_5k_R1.dex.fastq.gz
cat atac_v1_pbmc_5k_S1_L001_R3_001.dex.fastq.gz atac_v1_pbmc_5k_S1_L002_R3_001.dex.fastq.gz > atac_v1_pbmc_5k_R3.dex.fastq.gz

# index reference genome (bwa)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gzip -d hg19.fa.gz
snaptools index-genome --input-fasta=hg19.fa --output-prefix=hg19 --aligner=bwa --path-to-aligner=/usr/bin --num-threads=5

# Align (bwa)
snaptools align-paired-end \
--input-reference=../../Combine/3_Align/hg19.fa \
--input-fastq1=atac_v1_pbmc_5k_R1.dex.fastq.gz \
--input-fastq2=atac_v1_pbmc_5k_R3.dex.fastq.gz \
--output-bam=atac_v1_pbmc_5k.bam \
--aligner=bwa \
--read-fastq-command=zcat \
--min-cov=0 \
--num-threads=5 \
--if-sort=True \
--tmp-folder=./ \
--overwrite=TRUE




```

### Clustering Analysis

- R


```







```







