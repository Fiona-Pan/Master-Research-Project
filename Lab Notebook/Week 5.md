# WEEK 5

Zoom Meeting (2021.03.23 Tue):

- annotate peaks to genes/features
- extract top peaks based on ranked signal value
- start process cardiomyocyte analysis

## simple scRNA-seq clustering analysis

- Seurat


```
rm(list=ls())
gc()
# Seurat: scRNA-seq 5k PBMC
# scATAC-seq data: https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_pbmc_5k
# scRNA-seq data: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3_nextgem
# input:
#   - feature-by-cell mtx: barcodes.tsv.gz & features.tsv.gz & matrix.mtx.gz

#-----------------
# Dependency
#-----------------
library(dplyr)
library(Seurat)
library(patchwork)
#-----------------
# Read in scRNA-seq data
#-----------------
pbmc.data <- Read10X(data.dir ="/Users/apple/Downloads/Master_Research/data/5k_PBMC/5k_pbmc_v3_nextgem_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc5k", min.cells = 3, min.features = 100)
pbmc 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 30)

#-----------------
# Normalization
#-----------------
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))

#-----------------
# linear transformation & dimensional reduction
#-----------------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# DimPlot(pbmc, reduction = "pca")
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:20)

# ElbowPlot(pbmc)

#-----------------
# Clustering
#-----------------
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# saveRDS(pbmc, file = "5k_scRNA-seq_cluste.rds")

#-----------------
# Finding differentially expressed features
#-----------------
# find markers that define clusters via differential expression
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
#                                "CD8A"))
# 
# top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#-----------------
# Assigning cell type identity to clusters
#-----------------
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet","Monocytes","gdT","Neutrophil","T-reg")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



pbmc$celltype=Idents(pbmc)
ctypes <- as.vector(pbmc$celltype)
names(ctypes) <- names(pbmc$celltype)
pbmc <- AddMetaData(pbmc, metadata = ctypes, col.name = 'celltype')


saveRDS(pbmc, file = "pbmc5k_scRNA.rds")

```

![1](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/2.%20scRNA-seq-clustering/pdfs-1.jpg)
![2](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/2.%20scRNA-seq-clustering/pdfs-2.jpg)
![3](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/2.%20scRNA-seq-clustering/pdfs-3.jpg)



## scATAC-seq & scRNA-seq integration

```
# Integration of scATAC-seq & scRNA-seq
# annotate the single cell ATAC-seq clusters based on corresponding scRNA-seq dataset
# https://github.com/r3fang/SnapATAC/blob/master/examples/10X_PBMC_15K/README.md#annotation
# https://github.com/r3fang/SnapATAC/tree/master/examples/10X_brain_5k#gene_tsne
# scATAC-seq data: https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_v1_pbmc_5k
# scRNA-seq data: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3_nextgem

# input:
#   - scATAC profile: atac_v1_pbmc_5k.snap.rds
#   - scRNA profile: pbmc5k_scRNA.rds
#   - hg19 annotation: gencode.v19.annotation.gene.bed

#-----------------
# Dependency
#-----------------
library(Seurat)
library(GenomicRanges)
library(SnapATAC)
#-----------------
# Read in scATAC-seq and scRNA-seq data
#-----------------
x.sp = readRDS("atac_v1_pbmc_5k.snap.rds")
pbmc.rna=readRDS("pbmc5k_scRNA.rds")
pbmc.rna$tech = "rna";
variable.genes = VariableFeatures(object = pbmc.rna); #2000

genes.df = read.table("/Users/apple/Downloads/Master_Research/data/5k_PBMC/gencode.v19.annotation.gene.bed");
genes.gr = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4]);
genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)]; #1806

#-----------------
# Add cell-by-gene mtx
#-----------------
# reload the bmat
x.sp = addBmatToSnap(x.sp)

# Add cell x gene count matrix from scRNA-seq dataset
x.sp = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=3
)
x.sp #number of genes: 1806

#-----------------
# Seurat::transfer celltype labels
#-----------------
# convert the snap object to Seurat object 
#  identify anchors between the scATAC-seq dataset and the scRNA-seq dataset 
# and use these anchors to transfer the celltype labels 
pbmc.atac <- snapToSeurat(
  obj=x.sp, 
  eigs.dims=1:20, 
  norm=TRUE,
  scale=TRUE
)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc.rna, 
  query = pbmc.atac, 
  features = variable.genes, 
  reference.assay = "RNA", 
  query.assay = "ACTIVITY", 
  reduction = "cca"
)

# transfer the cluster ids, outputs a matrix with predictions and confidence scores
celltype.predictions <- TransferData(
  anchorset = transfer.anchors, 
  refdata = pbmc.rna$celltype,
  weight.reduction = pbmc.atac[["SnapATAC"]],
  dims = 1:20
);

x.sp@metaData$predicted.id = celltype.predictions$predicted.id;
x.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max);
x.sp@cluster = as.factor(x.sp@metaData$predicted.id);

#-----------------
# Create psudo multiomics cells
#-----------------
refdata <- GetAssayData(
  object = pbmc.rna, 
  assay = "RNA", 
  slot = "data"
)

imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = pbmc.atac[["SnapATAC"]], 
  dims = 1:20
)
x.sp@gmat = t(imputation@data)

#-----------------
#  Remove cells of low prediction score
#-----------------
hist(
  x.sp@metaData$predict.max.score, 
  xlab="prediction score", 
  col="lightblue", 
  xlim=c(0, 1),
  main="PBMC 10X"
)
abline(v=0.5, col="red", lwd=2, lty=2)
table(x.sp@metaData$predict.max.score > 0.5);
x.sp = x.sp[x.sp@metaData$predict.max.score > 0.5,];
x.sp

png("umap-cell-annotation-scATAC.png")
pdf("umap-cell-annotation-scATAC.pdf")
plotViz(
  obj=x.sp,
  method="umap", 
  main="PBMC 10X",
  point.color=x.sp@metaData[,"predicted.id"], 
  point.size=0.5, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)
dev.off()
#-----------------
# Gene expression projected onto UMAP
#-----------------

marker.genes = c(
  "IL32", "LTB", "CD3D",
  "IL7R", "LDHB", "FCGR3A", 
  "CD68", "MS4A1", "GNLY", 
  "CD3E", "CD14", "CD14", 
  "FCGR3A", "LYZ", "PPBP", 
  "CD8A", "PPBP", "CST3", 
  "NKG7", "MS4A7", "MS4A1", 
  "CD8A"
)
par(mfrow = c(3, 3))
for(i in 1:9){
  j = which(colnames(x.sp@gmat) == marker.genes[i])
  
  plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@gmat[,j],
    method="umap", 
    main=marker.genes[i],
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0.01, 0.99))
  }
dev.off()
#-----------------
# Identify differentially accessible peaks
#-----------------

clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)]

idy.ls = lapply(clusters.sel, function(cluster_i){
  DARs = findDAR(
    obj=x.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    bcv=0.4,
    test.method="exactTest",
    seed.use=10
  );
  DARs$FDR = p.adjust(DARs$PValue, method="BH");
  idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  if((x=length(idy)) < 2000L){
    PValues = DARs$PValue;
    PValues[DARs$logFC < 0] = 1;
    idy = order(PValues, decreasing=FALSE)[1:2000];
    rm(PValues); # free memory
  }
  idy
})
names(idy.ls) = clusters.sel;

pdf("differential-accessible-peaks.pdf")
par(mfrow = c(3, 3));
for(cluster_i in clusters.sel){
  print(cluster_i)
  idy = idy.ls[[cluster_i]];
  vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
  vals.zscore = (vals - mean(vals)) / sd(vals);
  plotFeatureSingle(
    obj=x.sp,
    feature.value=vals.zscore,
    method="tsne", 
    main=cluster_i,
    point.size=0.1, 
    point.shape=19, 
    down.sample=5000,
    quantiles=c(0.01, 0.99)
  );
}
dev.off()

saveRDS(x.sp, file="atac_v1_pbmc_5k.final.rds")

```

![4](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/3.%20scATAC-scRNA-integration/pdfs-1.jpg)
![5](https://github.com/Fiona-Pan/Master-Research-Project/blob/main/plots/3.%20scATAC-scRNA-integration/pdfs-2.jpg)


