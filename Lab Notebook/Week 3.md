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

# pre-processing
snaptools snap-pre \
--input-file=atac_v1_pbmc_5k.bam \
--output-snap=atac_v1_pbmc_5k.snap \
--genome-name=hg19 \
--genome-size=hg19.chrom.sizes \
--min-mapq=30 \
--min-flen=0 \
--max-flen=1000 \
--keep-chrm=TRUE \
--keep-single=FALSE \
--keep-secondary=FALSE \
--overwrite=True \
--min-cov=100 \
--verbose=True

# Cell-by-bin matrix
snaptools snap-add-bmat \
--snap-file=atac_v1_pbmc_5k.snap \
--bin-size-list 1000 2000 5000 10000 \
--verbose=True


```

### Clustering Analysis

- R


```







```







