# WEEK 3 & WEEK 4


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
# input: 
#   - snap file(converted from BAM by snapTools): atac_v1_pbmc_5k.snap
#   - per barcode metrics: atac_v1_pbmc_5k_singlecell.csv
#   - hg19 blacklist: Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz
rm(list=ls())
gc()
#-----------------
# Dependencies
#-----------------
library(SnapATAC)
library(viridisLite)
library(ggplot2)
library(GenomicRanges)
library(leiden)
library(umap)

#-----------------
# Barcode Selection
#-----------------
# select high-quality barcodes based on 
  # 1) # of unique fragments; 
  # 2) fragments in promoter ratio
x.sp=createSnap(
  file="/Users/apple/Downloads/Master_Research/data/5k_PBMC/atac_v1_pbmc_5k.snap",
  sample="atac_v1_pbmc_5k",
  num.cores=1
)
x.sp #number of barcodes: 91652
summarySnap(x.sp)

plotBarcode(
  obj=x.sp, 
  pdf.file.name=NULL, 
  pdf.width=7, 
  pdf.height=7, 
  col="grey",
  border="grey",
  breaks=50
)

barcodes=read.csv("/Users/apple/Downloads/Master_Research/data/5k_PBMC/atac_v1_pbmc_5k_singlecell.csv",
                  head=TRUE)
barcodes$barcode=substr(barcodes$barcode,1,16)
barcodes=barcodes[2:nrow(barcodes),]

promoter_ratio=(barcodes$promoter_region_fragments+1)/(barcodes$passed_filters+1)
UMI=log(barcodes$passed_filters+1,10) #unique molecular identifier 
data=data.frame(UMI=UMI,promoter_ratio=promoter_ratio)
barcodes$promoter_ratio=promoter_ratio
barcodes$logUMI=UMI

#plot
p1=ggplot(
  data,aes(x=UMI,y=promoter_ratio)) +
  geom_point(size=0.1,col="grey") +
  theme_classic() +
  ggtitle("5X PBMC") +
  ylim(0,1) + xlim(0,6) +
  labs(x="log10(UMI)",y="promoter ratio")
p1

# filter
barcodes.sel=barcodes[which(UMI>=3 & UMI <=5 & promoter_ratio >=0.15 & promoter_ratio <=0.8),]
rownames(barcodes.sel)=barcodes.sel$barcode
x.sp=x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),]
x.sp@metaData=barcodes.sel[x.sp@barcode,]
x.sp # number of barcodes: 4841 (was 91652)

```

UMI: unique molecular identifier 

![1](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/Rplot-barcodes%20filtration(Page1).jpg)

```
#-----------------
# Add cell-by-bin matrix
#-----------------
# add the cell-by-bin matrix to snap object
showBinSizes("/Users/apple/Downloads/Master_Research/data/5k_PBMC/atac_v1_pbmc_5k.snap") #1000  2000  5000 10000
x.sp=addBmatToSnap(x.sp,bin.size=5000,num.cores=1) #5kb resolution
x.sp #number of barcodes: 4841; number of bins: 627478

#-----------------
# Matrix binarization
#-----------------
# convert cell-by-bin coount matrix to a binary matrix
# remove 0.1% items of the highest coverage in the count matrix, convert remaining non-zero to 1
x.sp=makeBinary(x.sp,mat="bmat")

#-----------------
# Bin Filtering
#-----------------
# filter out bins overlapping with the ENCODE blacklist
# prevent potential artifacts
blacklist=read.table("/Users/apple/Downloads/Master_Research/data/5k_PBMC/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz")
blacklist.gr=GRanges(
  blacklist[,1],
  IRanges(blacklist[,2],blacklist[,3])
)
idy=queryHits(findOverlaps(x.sp@feature,blacklist.gr))
if(length(idy)>0){x.sp=x.sp[,-idy,mat="bmat"]}
x.sp #number of bins: 625212 (was: 627478)

# remove unwanted chromosomes
chr.exclude=seqlevels(x.sp@feature)[grep("random|chrM",seqlevels(x.sp@feature))]
idy=grep(paste(chr.exclude,collapse="|"),x.sp@feature)
if(length(idy)>0){x.sp=x.sp[,-idy,mat="bmat"]}
x.sp #number of bins: 624711 (was: 625212)

# remove top 5% bins that overlap with invariant features (ex. promoters of house keeping genes)
bin.cov=log10(Matrix::colSums(x.sp@bmat)+1)
hist(
  bin.cov[bin.cov>0],
  xlab="log10(bin cov)",
  main="log10(Bin Cov)",
  col="lightblue",
  xlim=c(0,5)
)
bin.cutoff=quantile(bin.cov[bin.cov>0],0.95)
idy=which(bin.cov <= bin.cutoff & bin.cov >0)
x.sp=x.sp[,idy,mat="bmat"]
x.sp #number of bins: 531640 (was: 624711)

#  further remove any cells of bin coverage less than 1,000
idx = which(Matrix::rowSums(x.sp@bmat) > 1000)
x.sp = x.sp[idx,]
x.sp #number of barcodes: 4366 (was: 4841)

#-----------------
# Dimensionality reduction
#-----------------
# compute diffusion maps for dimentionality reduction
row.covs.dens <- density(
  x = x.sp@metaData[,"logUMI"], 
  bw = 'nrd', adjust = 1
)
sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps)
set.seed(1)
idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = 4000, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
x.query.sp = x.sp[-idx.landmark.ds,];

x.landmark.sp = runDiffusionMaps(
  obj= x.landmark.sp,
  input.mat="bmat", 
  num.eigs=50
);
x.landmark.sp@metaData$landmark = 1;
x.query.sp = runDiffusionMapsExtension(
  obj1=x.landmark.sp, 
  obj2=x.query.sp,
  input.mat="bmat"
);
x.query.sp@metaData$landmark = 0;
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp@metaData["sample"] = x.sp@sample
x.sp = x.sp[order(x.sp@metaData[,"sample"])];

#-----------------
# Determine significant components
#-----------------
plotDimReductElbow(
  obj=x.sp, 
  point.size=1.5,
  point.shape=19,
  point.color="red",
  point.alpha=1,
  pdf.file.name=NULL,
  pdf.height=7,
  pdf.width=7,
  labs.title="PCA Elbow plot",
  labs.subtitle=NULL
)

plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
)


df_out = data.frame(t(x.sp@smat@dmat)[1:20,],row.names = sprintf("PC_%s",seq(1:20)))
colnames(df_out) = as.character(x.sp@metaData$barcode)
dim(df_out)
df_out[1:5,1:5]
```

![2](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/Rplot-PCA%20dimension%20selection-Elbow%20point(Page1).jpg)


![3](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/Rplot-PCA%20dimension%20selection(Page1).jpg)


```
#-----------------
# Graph-based clustering
#-----------------
# construct K Nearest Neighbor Graph
x.sp=runKNN(
  obj=x.sp,
  eigs.dims=1:20, #20 dims
  k=15
)

x.sp=runCluster(
  obj=x.sp,
  tmp.folder = tempdir(),
  louvain.lib="leiden", #"R-igraph"
  seed.use=10
)



#x.sp@metaData$cluster=x.sp@cluster
#-----------------
# Visualization
#-----------------
# umap
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="umap",
  seed.use=10
);

par(mfrow = c(2, 2))
plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@cluster, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",,
  down.sample=10000,
  legend.add=FALSE
);

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData[,"logUMI"],
  method="umap", 
  main="Read Depth",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotViz(
  obj= x.sp,
  method="umap", 
  main="Sample",
  point.size=0.2, 
  point.shape=19, 
  point.color=x.sp@sample, 
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
);

plotViz(
  obj= x.sp,
  method="umap", 
  main="Landmark",
  point.size=0.2, 
  point.shape=19, 
  point.color=x.sp@metaData[,"landmark"], 
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
);

#t-SNE
x.sp=runViz(
  obj=x.sp,
  tmp.folder = tempdir(),
  dims=2,
  eigs.dims=1:20,
  method="Rtsne",
  seed.use=10
)

par(mfrow = c(2, 2))
plotViz(
  obj=x.sp,
  method="tsne",
  main="5K PBMC",
  point.color=x.sp@cluster,
  point.size=1,
  point.shape=19,
  point.alpha=0.8,
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)

plotFeatureSingle(
  obj=x.sp,
  feature.value = log(x.sp@metaData[,"passed_filters"]+1,10),
  method="tsne",
  main="5K PBMC read depth",
  point.size=0.2,
  point.shape=19,
  down.sample=10000,
  quantiles=c(0.01,0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters,
  method="tsne", 
  main="5K PBMC FRiP",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99) # remove outliers
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$duplicate / x.sp@metaData$total,
  method="tsne", 
  main="5K PBMC Duplicate",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99) # remove outliers
)

#-----------------
# Gene based annotation
#-----------------
# annotate identified cell clusters
# create cell-by-gene matrix and visualize the enrichment of marker genes


#-----------------
# Heretical clustering
#-----------------
# cells belong to the same cluster are pooled to create the aggregate signal
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})

hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");

plotViz(
  obj=x.sp,
  method="tsne", 
  main="5K PBMC Cluster",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)

plot(hc,hang=-1,xlab="")

```

![4](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/Rplot-UMAP%2BtSNE%20clustering(Page1).jpg)

![5](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/Rplot-tSNE%20clustering(Page1).jpg)

![6](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/Rplot-UMAP%20clustering(Page1).jpg)


```
#-----------------
# Identify peaks
#-----------------
# aggregate cells from each cluster to create an ensemble track for peak calling
# snaptools: /Users/apple/miniconda3/bin/snaptools
# macs2: /Users/apple/miniconda3/bin/macs2
# runMACS(
#   obj=x.sp[which(x.sp@cluster==1),], #first cluster
#   output.prefix="atac_v1_pbmc_5k.1",
#   path.to.snaptools = "/Users/apple/miniconda3/bin/snaptools",
#   path.to.macs = "/Users/apple/miniconda3/bin/macs2",
#   gsize="hs",
#   buffer.size=500,
#   num.cores = 5,
#   macs.options = "--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
#   tmp.folder = tempdir()
# )
# 
# clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)]
# # for all clusters with more than 150 cells
# peaks.ls = mclapply(seq(clusters.sel), function(i){
#   print(clusters.sel[i]);
#   runMACS(
#     obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
#     output.prefix=paste0("atac_v1_pbmc_5k.", gsub(" ", "_", clusters.sel)[i]),
#     path.to.snaptools="/Users/apple/miniconda3/bin/snaptools",
#     path.to.macs="/Users/apple/miniconda3/bin/macs2",
#     gsize="hs",
#     buffer.size=500, 
#     num.cores=1,
#     macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
#     tmp.folder=tempdir()
#   );
# }, mc.cores=5)

peaks.names = system("ls | grep narrowPeak", intern=TRUE);

peak.gr.ls = lapply(peaks.names, function(x){ # a list of GRanges of clusters
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})

peak.gr = reduce(Reduce(c, peak.gr.ls)) # reduce into one GRange object for all clusters
peak.gr


#-----------------
# Create cell-by-peak matrix
#-----------------
# 
peaks.df = as.data.frame(peak.gr)[,1:3]
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

saveRDS(x.sp, file="atac_v1_pbmc_5k.snap.rds")
```


![7](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/IGV-plot-by-clusters.png)

create cell-by-peak matrix and add to the snap file
```
# 
  # snaptools snap-add-pmat \
    # --snap-file /Users/apple/Downloads/Master_Research/data/5k_PBMC/atac_v1_pbmc_5k.snap \
    # --peak-file /Users/apple/Documents/peaks.combined.bed 


```



```
#-----------------
# Add cell-by-peak matrix
#-----------------
#x.sp = readRDS("atac_v1_pbmc_5k.snap.rds")
x.sp = addPmatToSnap(x.sp)
x.sp = makeBinary(x.sp, mat="pmat")
x.sp #number of peaks: 102613

#-----------------
table(x.sp@cluster) # number of cells in each cluster
x.sp@pmat #cell-by-peak matrix
x.sp@pmat@Dim #4366 (cells) 102613 (peaks) 
length(x.sp@pmat@x)/length(x.sp@pmat) # 0.07368226 (7.4% of non-zero count peaks)

# pmat cluster1
cluster1=x.sp[which(x.sp@cluster==1),]
cluster1@pmat
cluster1@pmat@Dim #571 102613
length(cluster1@pmat@x)/length(cluster1@pmat) #0.06124419

# pmat by-cluster
clusters.sel = names(table(x.sp@cluster))
x.sp.cluster.ls=list()
x.sp.cluster.names=c()
for(i in seq(clusters.sel)){
  names=sprintf("cluster%s",i)
  x.sp.clusters=x.sp[which(x.sp@cluster==i),]
  x.sp.cluster.ls=append(x.sp.cluster.ls,x.sp.clusters)
  x.sp.cluster.names=append(x.sp.cluster.names,names)
}
names(x.sp.cluster.ls)<-x.sp.cluster.names
x.sp.cluster.ls

for(i in 1:length(x.sp.cluster.ls)){
  cluster_i=x.sp.cluster.ls[[i]]
  names=names(x.sp.cluster.ls)[i]
  nz_peaks=length(cluster_i@pmat@x)
  total_length=length(cluster_i@pmat)
  nz_counts=nz_peaks/total_length
  #print(sprintf("%s: %s cells & %s peaks",names,cluster_i@pmat@Dim[1],cluster_i@pmat@Dim[2]))
  print(sprintf("non-zero peaks couts: %s: %s",names,round(nz_counts,4)))
  }


sum(x.sp.cluster.ls$cluster1@pmat[1,]) # No. peaks in one cell (row)
sum(x.sp.cluster.ls$cluster1@pmat[,1]) # No. cells in one peak (col)

pk_in_cell_cluster_1=c() # how many peaks in each cell
cell_in_pk_cluster_1=c() # how many cells in each peak

for(i in 1:nrow(x.sp.cluster.ls$cluster1@pmat)){
  cell_name=rownames(x.sp.cluster.ls$cluster1@pmat)[i]
  pk_in_cell=sum(x.sp.cluster.ls$cluster1@pmat[i,])
  pk_in_cell_cluster_1=append(pk_in_cell_cluster_1,sprintf("cell %s: %s peaks total",i,pk_in_cell))
  
  cell_in_pk=sum(x.sp.cluster.ls$cluster1@pmat[,i])
  cell_in_pk_cluster_1=append(cell_in_pk_cluster_1,sprintf("peak %s: %s cells",i,cell_in_pk))
  }

#-----------------

```

**Worth Mention:**

1. Cluster info

```
table(x.sp@cluster) 
number of cells in each cluster
  1   2   3   4   5   6   7   8   9  10  11  12  13  14 
663 611 505 489 398 368 241 208 193 172 160 134 122 102 
```

2. Peak info

```
x.sp@pmat # peak-by-cell matrix
x.sp@pmat@Dim # 4366 (cells) 102613 (peaks) 
length(x.sp@pmat@x)/length(x.sp@pmat) # 0.07368226 (7.4% of non-zero count peaks)

by-cluster: # non-zero count peaks
[1] "non-zero peaks couts: cluster1: 0.0619"
[1] "non-zero peaks couts: cluster2: 0.0596"
[1] "non-zero peaks couts: cluster3: 0.1007"
[1] "non-zero peaks couts: cluster4: 0.0826"
[1] "non-zero peaks couts: cluster5: 0.0829"
[1] "non-zero peaks couts: cluster6: 0.0658"
[1] "non-zero peaks couts: cluster7: 0.0807"
[1] "non-zero peaks couts: cluster8: 0.0725"
[1] "non-zero peaks couts: cluster9: 0.0749"
[1] "non-zero peaks couts: cluster10: 0.0607"
[1] "non-zero peaks couts: cluster11: 0.0687"
[1] "non-zero peaks couts: cluster12: 0.0787"
[1] "non-zero peaks couts: cluster13: 0.0699"
[1] "non-zero peaks couts: cluster14: 0.061"
```



