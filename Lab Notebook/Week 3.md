# WEEK 3


Zoom meeting (2020.03.09 Tue):

- read pre-processed dataset in R (peak-by-cell matrix)
- clustering and find biological meanings

### Snaptools Pipeline

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
--output-fastq=atac_v1_pbmc_5k_S1_L002_R3_002.dex.fastq.gz \
--index-fastq-list atac_v1_pbmc_5k_S1_L002_R2_001.fastq.gz

# combine lanes
cat atac_v1_pbmc_5k_S1_L001_R1_001.dex.fastq.gz atac_v1_pbmc_5k_S1_L002_R1_001.dex.fastq.gz > atac_v1_pbmc_5k_R1.dex.fastq.gz
cat atac_v1_pbmc_5k_S1_L001_R3_001.dex.fastq.gz atac_v1_pbmc_5k_S1_L002_R3_001.dex.fastq.gz > atac_v1_pbmc_5k_R3.dex.fastq.gz

# index reference genome (bwa)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gzip -d hg19.fa.gz
snaptools index-genome --input-fasta=hg19.fa --output-prefix=hg19 --aligner=bwa --path-to-aligner=/usr/bin --num-threads=5

# Align (bwa)





```

