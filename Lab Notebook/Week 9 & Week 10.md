# WEEK 9 & WEEK 10

##### Table of Contents  
[Headers](#headers)  
[Emphasis](#emphasis)  
   
   
<a name="headers"/>
## Headers
xxx

***


**Zoom Meeting:**

- downstream scATAC-seq analysis:
  - integration of scRNA-seq dataset
  - peak annotation, gene signatures
  - enrichment analysis (e.g enriched GO term for each cluster)
  - investigation of cell population changes b/t control and kD
  - ...

- sample ID
  - Library_20190731_001: Control sample (sometimes referred to as 0XAV)
  - Library_20190731_002: TMEM88 KD sample (sometimes referred to as Dox)

**Previous**

- clustered by bins (UMAP + tSNE)


### 0. clustered by peaks (UMAP)

- input
  - peak-by-cell matrix


- output
  - L1_pmat_clustering.snap.RDS (SnapATAC Obj.)
  - cluster by peaks (UMAP)

```
#-----------------
# L1: cluster by peaks & save RDS
#-----------------
library(leiden)
library(umap)
setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS")

x.sp<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/Library_20190731_001_S2_after_add_peaks.snap.rds")
x.sp
row.covs.dens <- density(
  x = x.sp@metaData[,"logUMI"], 
  bw = 'nrd', adjust = 1
);
sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps);
idx.landmark.ds <- sort(sample(x = seq(nrow(x.sp)), size = 6000, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
x.query.sp = x.sp[-idx.landmark.ds,];

x.landmark.sp = runDiffusionMaps(
  obj= x.landmark.sp,
  input.mat="pmat", 
  num.eigs=50
)
x.landmark.sp@metaData$landmark = 1
#------!!!-----------
x.query.sp = runDiffusionMapsExtension(
  obj1=x.landmark.sp, 
  obj2=x.query.sp,
  input.mat="pmat"
)
#-------!!!----------
x.query.sp@metaData$landmark = 0;
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp = x.sp[order(x.sp@metaData[,"sample"])]; 
x.sp

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

x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
)

x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  seed.use=10,
  resolution=0.7
)

x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="umap",
  seed.use=10
)
pdf("L1-pmat-cluster.pdf")
plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@cluster, 
  point.size=0.2, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)
dev.off()

saveRDS(x.sp,"L1_pmat_clustering.snap.RDS")
```

### 1. scRNA-seq operation

- input: 
  - raw scRNA-seq metadata
  - raw scRNA-seq expression matrix
  - three samples (need control (0XAV) and TMEM88 KD sample (Dox))

- output:
  - scrna1-Seurat.rds (Seurat Obj.)

```
rm(list=ls())
gc()
#-----------------
# Dependencies
#-----------------
setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/scRNA")

library(dplyr)
library(Seurat)
library(patchwork)
library(crowbar)
library(Matrix)
library(plyr)
#-----------------
# Read in files
#-----------------
# raw metadata and raw exp contain three samples
# only need control (0XAV) and TMEM88 KD sample (Dox)
rawmetadata<-readRDS("metadata.RDS")
rawmetadata$cluster<-NULL
rawexp<-readRDS("gmat_exp_mtx.RDS")

metadata<-rawmetadata[rawmetadata$sample==c("0Xav","Dox"),]
exp<-rawexp[rawexp$sample==c("0Xav","Dox"),]

metadata1<-rawmetadata[rawmetadata$sample=="0Xav",]


exp1<-rawexp[rawexp$sample=="0Xav",]


gmat1<-exp1[,8:ncol(exp1)]
gmat1<-as.data.frame(gmat1)
rownames(gmat1)<-exp1$cell_barcode
gmat1<-as.matrix(gmat1,sparse=TRUE)
gmat1<-as(gmat1,"dgCMatrix")
gmat1<-t(gmat1)

#-----------------
# create Seurat Object
#-----------------
scrna1<-CreateSeuratObject(counts = gmat1, project = "scRNA", min.cells = 3, min.features = 100)

scrna1[["seurat_clusters"]]<-metadata1$clust0.2

# 0    1    2    3    4    5    6    7    8    9   10 
# 3459 3753  868 1449 2909  877 1121  592  537  211  233 

# 0=Lateral Plate mesoderm 
# 1=Definitive Endoderm
# 2=Mesoderm
# 3=Cardiac cells
# 4=Posterior foregut endoderm
# 5=Mesendoderm
# 6=Anterior foregut endoderm
# 7=Paraxial mesoderm
# 8=Endocardial Endothelium
# 9=Axial Mesoderm
# 10=MT1F+ cells


celltype<-revalue(metadata1$clust0.2, c("0"="Lateral Plate mesoderm", 
                    "1"="Definitive Endoderm",
                    "2"="Mesoderm",
                    "3"="Cardiac cells",
                    "4"="Posterior foregut endoderm",
                    "5"="Mesendoderm",
                    "6"="Anterior foregut endoderm",
                    "7"="Paraxial mesoderm",
                    "8"="Endocardial Endothelium",
                    "9"="Axial Mesoderm",
                    "10"="MT1F+ cells"))


scrna1[["celltype"]]<-celltype
scrna1@active.ident<-scrna1$celltype

scrna1<-FindVariableFeatures(scrna1, selection.method = "vst", nfeatures = 2000)

umap<-data.frame(UMAP_1=metadata1$UMAP_1,UMAP_2=metadata1$UMAP_2)
rownames(umap)<-metadata1$cell_barcode

scrna1[["UMAP_1"]]<-umap$UMAP_1
scrna1[["UMAP_2"]]<-umap$UMAP_2


scrna1@assays$RNA@scale.data<-as.matrix(gmat1)
scrna1 <- RunPCA(scrna1, features = VariableFeatures(object = scrna1))
scrna1<-RunUMAP(scrna1, dims = 1:8)

# all.genes <- rownames(scrna1)
# scrna1 <- ScaleData(scrna1, features = all.genes)
# scrna1 <- RunPCA(scrna1, features = VariableFeatures(object = scrna1))
# scrna1<-RunUMAP(scrna1, dims = 1:10)

scrna1@reductions$umap@cell.embeddings<-as.matrix(umap)

scrna1.markers <- FindAllMarkers(scrna1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scrna1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- scrna1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#-----------------
# Plots
#-----------------
# VlnPlot(pbmc.rna, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# 
# plot1 <- VariableFeaturePlot(scrna1)
# plot2 <- LabelPoints(plot = plot1, points = top10$gene, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))
# plot2

# VizDimLoadings(scrna1, dims = 1:2, reduction = "pca")

# pdf("1-UMAP-metadata.pdf")
DimPlot(scrna1, reduction = "umap")
DimPlot(scrna1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# dev.off()
#-----------------
# Results
#-----------------

# Highly Variable Features
# HVFInfo(scrna1)
# top10 <- head(VariableFeatures(scrna1), 10)

#-----------------
# Save
#-----------------
saveRDS(scrna1, file = "scrna1-Seurat.rds")
```

### 2. scATAC-seq & scRNA-seq integration

- input:
  - scATAC-seq: peak-by-cell matrix
  - scRNA-seq: Seurat Obj. 
  - metadata: Library_20190731_001_singlecell.csv, Library_20190731_001_fragments.tsv.gz

- output:
  - L1-scATAC-integration-Seurat.rds (Seurat Obj.)
  - integration clustering (UMAP)

```

#-----------------
# Signac(Seurat): peak anno
#-----------------
rm(list=ls())
#-----------------
# Dependencies
#-----------------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(stringr)

setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte")
#-----------------
# File Prep
#-----------------
x.sp<-readRDS("RDS/Library_20190731_001_S2_after_add_peaks.snap.rds")

# pmat
pmat<-x.sp@pmat
peak_names<-x.sp@peak@elementMetadata$name
a<-str_replace(peak_names,"b'","")
peak_names<-str_replace(a,"'","")
colnames(pmat)<-peak_names
rownames(pmat)<-paste0(rownames(pmat),"-1")
pmat_t<-t(pmat)
counts<-pmat_t

# metadata
metadata <- read.csv(
  file = "Library_20190731_001_singlecell.csv",
  header = TRUE,
  row.names = 1
)

#-----------------
# Seurat Obj.
#-----------------
# creating a Seurat object using the peak/cell matrix and cell metadata generated by cellranger-atac
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'Library_20190731_001_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
atac #279423 features across 7010 samples within 1 assay
atac[["peaks"]]
granges(atac)




annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
Annotation(atac) <- annotations

# atac <- NucleosomeSignal(object = atac)
# atac <- TSSEnrichment(object = atac, fast = FALSE)
# 
# atac$high.tss <- ifelse(atac$TSS.enrichment > 2, 'High', 'Low')
# TSSPlot(atac, group.by = 'high.tss') + NoLegend()
# 
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac)
# DepthCor(atac)

atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 1:20)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 1:20)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)
DimPlot(object = atac, label = TRUE) + NoLegend()
cluster<-Idents(atac)
table(cluster)

# 0   1   2   3   4   5   6   7   8   9  10  11  12  13 
# 916 875 701 671 623 601 517 403 395 324 285 277 266 156

gene.activities <- GeneActivity(atac)
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)

DefaultAssay(atac) <- 'RNA'

rna<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/scRNA/scrna1-Seurat.rds")
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$celltype,
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

atac <- AddMetaData(object = atac, metadata = predicted.labels)
atac$predicted.id[atac$predicted.id=="MT1F+ cells"]<-"MT1F+"
atac$predicted.id[atac$predicted.id=="Anterior foregut endoderm"]<-"An.for.en"
atac$predicted.id[atac$predicted.id=="Axial Mesoderm"]<-"Ax.Me"
atac$predicted.id[atac$predicted.id=="Cardiac cells"]<-"Car.cell"
atac$predicted.id[atac$predicted.id=="Definitive Endoderm"]<-"Def.En"
atac$predicted.id[atac$predicted.id=="Endocardial Endothelium"]<-"Endo.Endo"
atac$predicted.id[atac$predicted.id=="Lateral Plate mesoderm"]<-"Lat.P.me"
atac$predicted.id[atac$predicted.id=="Mesendoderm"]<-"Mesen"
atac$predicted.id[atac$predicted.id=="Paraxial mesoderm"]<-"Pa.meso"
atac$predicted.id[atac$predicted.id=="Posterior foregut endoderm"]<-"Pps.for.endo"
atac$predicted.id[atac$predicted.id=="Mesoderm"]<-"Meso"


plot1 <- DimPlot(
  object = rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot1.1 <- DimPlot(
  object = rna,
  group.by = 'celltype',
  label = FALSE,
  repel = TRUE)  + ggtitle('scRNA-seq')


plot2 <- DimPlot(
  object = atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle('scATAC-seq')

plot2.1 <- DimPlot(
  object = atac,
  group.by = 'predicted.id',
  label = FALSE,
  repel = TRUE) + ggtitle('scATAC-seq')

plot3 <- DimPlot(
  object = atac,
  label = FALSE,
  repel = TRUE) + ggtitle('scATAC-seq')

plot3.1 <- DimPlot(
  object = atac,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pdf("L1-rna-atac.pdf")
plot1 + plot2
dev.off()

pdf("L1-rna-atac-v2.pdf")
plot1.1 + plot2
dev.off()

pdf("L1-atac-atac.pdf")
plot2 + plot3
dev.off()

pdf("L1-atac-atac-v2.pdf")
plot2.1 + plot3.1
dev.off()


saveRDS(atac,"/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/L1-scATAC-integration-Seurat.rds")

```
