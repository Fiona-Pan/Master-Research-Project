# Week 11

**Zoom Meeting**

- merge peak-by-cell matrix of the scATAC-seq control and KD samples together 
- ideas of find out genes/peaks/biological differences from for scATAC cluster that were predicted as one large cell type group from scRNA integration
  - Sophie suggest to investigate clusters that were predicted to be cardiac (maybe  lateral plate mesoderm in the KD condition)
- **Find differentially accessible peaks between clusters**
  - perform differential accessibility (DA) test utilizing logistic regression or look at fold change accessibility between two group of cells
  - find closest gene to each of differentially accessible peaks followed by gene ontology enrichment analysis on those gene sets
  - motif variability analysis & de novo motif discovery to identify master regulators that are enriched in DARs
  - Plotting genomic regions
- **sample-specific: pathway-enrichment analysis**
  - Hypothesis: elevated Wnt signalling in the *TMEM88* KD condition

# Merge scATAC

### **Seurat:: [Merging Objects](https://satijalab.org/signac/articles/merging.html)**

```
#-----------------#-----------------
# # Seurat:: Merging Objects
# # https://satijalab.org/signac/articles/merging.html
# # https://satijalab.org/seurat/archive/v3.0/merge_vignette.html
#   - 1. Merge from pmat
#   - 2. Merge from FeatureMatrix
#-----------------#-----------------
# rm(list=ls())
# #----------------
# # 1. Merge from pmat
# # "Merge-by-L1-L2-pmat-Seurat.Obj.rds"
# #-----------------
# library(Signac)
# library(Seurat)
# library(GenomicRanges)
# library(future)
# 
# p1<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/pmat/L1-pmat.rds")
# p2<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/pmat/L2-pmat.rds")
# 
# rownames(p1)<-gsub(pattern = '_', replacement = '', x = rownames(p1))
# rownames(p2)<-gsub(pattern = '_', replacement = '', x = rownames(p2))
# 
# 
# s1<-read.table(
#   file = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_001_singlecell.csv", #file/Library_20190731_001_singlecell.csv
#   stringsAsFactors = FALSE,
#   sep = ",",
#   header = TRUE,
#   row.names = 1
# )[-1, ]
# s1<-s1[s1$passed_filters > 100, ]
# 
# s2<-read.table(
#   file = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_002_singlecell.csv", #file/Library_20190731_002_singlecell.csv
#   stringsAsFactors = FALSE,
#   sep = ",",
#   header = TRUE,
#   row.names = 1
# )[-1, ]
# s2<-s2[s2$passed_filters > 100, ]
# 
# f1 <- CreateFragmentObject(
#   path = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_001_fragments.tsv.gz", #file/Library_20190731_001_fragments.tsv.gz
#   cells = rownames(s1)
# )
# 
# f2 <- CreateFragmentObject(
#   path = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_002_fragments.tsv.gz", #file/Library_20190731_002_fragments.tsv.gz
#   cells = rownames(s2)
# )
# 
# c1<-CreateChromatinAssay(counts = p1, sep = c(":", "-"),fragments = f1)
# counts1<-CreateSeuratObject(c1, assay = "ATAC",meta.data = s1)
# 
# 
# c2<-CreateChromatinAssay(counts = p2, sep = c(":", "-"),fragments = f2)
# counts2<-CreateSeuratObject(c2, assay = "ATAC",meta.data = s2)
# 
# 
# counts1$dataset <- 'WT'
# counts2$dataset <- 'KD'
# 
# 
# GetAssayData(counts1,slot="counts")[1:5,1:5]
# GetAssayData(counts2,slot="counts")[1:5,1:5]
# 
# counts1<-NormalizeData(counts1)
# counts2<-NormalizeData(counts2)
# 
# 
# combined <- merge(
#   x = counts1,
#   y = counts2,
#   add.cell.ids = c("WT","KD")
# )
# 
# 
# combined <- RunTFIDF(combined)
# combined <- FindTopFeatures(combined, min.cutoff = 20)
# combined <- RunSVD(combined)
# combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
# 
# pdf("merge-by-pmat.Seurat.pdf")
# DimPlot(combined)
# DimPlot(combined, group.by = 'dataset')
# dev.off()
# 
# saveRDS(combined,"Merge-by-L1-L2-pmat-Seurat.Obj.rds")

# #----------------
# # 2. Merge from FeatureMatrix
# # "Merge-by-FeatureMatrix.rds"
# #-----------------
# 
# library(Signac)
# library(Seurat)
# library(GenomicRanges)
# library(future)
# library(stringr)
# library(Matrix)
# 
# d1<-read.csv("../CSV/L1-peak-df.csv")
# d2<-read.csv("../CSV/L2-peak-df.csv")
# 
# g1<-GRanges(d1$chr,IRanges(d1$start,d1$end))
# g2<-GRanges(d2$chr,IRanges(d2$start,d2$end))
# 
# gc<-reduce(c(g1,g2))
# gcc<-as.data.frame(gc)
# pc<-str_c(gcc$seqnames,":",gcc$start,"-",gcc$end)
# 
# p1<-readRDS("../CSV/L1-pmat.rds")
# p2<-readRDS("../CSV/L2-pmat.rds")
# 
# p1<-p1[which(rownames(p1) %in% pc),]
# p2<-p2[which(rownames(p2) %in% pc),]
# 
# dim(p1)#[1] 131049   7010
# 
# dim(p2)#[1] 142551   4830
# 
# rownames(p1)<-gsub(pattern = '_', replacement = '', x = rownames(p1))
# rownames(p2)<-gsub(pattern = '_', replacement = '', x = rownames(p2))
# 
# s1<-read.table(
#   file = "../file/Library_20190731_001_singlecell.csv", #file/Library_20190731_001_singlecell.csv
#   stringsAsFactors = FALSE,
#   sep = ",",
#   header = TRUE,
#   row.names = 1
# )[-1, ]
# s1<-s1[s1$passed_filters > 100, ]
# 
# 
# s2<-read.table(
#   file = "../file/Library_20190731_002_singlecell.csv", #file/Library_20190731_002_singlecell.csv
#   stringsAsFactors = FALSE,
#   sep = ",",
#   header = TRUE,
#   row.names = 1
# )[-1, ]
# s2<-s2[s2$passed_filters > 100, ]
# 
# f1 <- CreateFragmentObject(
#   path = "../file/Library_20190731_001_fragments.tsv.gz", #file/Library_20190731_001_fragments.tsv.gz
#   cells = rownames(s1)
# )
# 
# f2 <- CreateFragmentObject(
#   path = "../file/Library_20190731_002_fragments.tsv.gz", #file/Library_20190731_002_fragments.tsv.gz
#   cells = rownames(s2)
# )
# 
# 
# p1<-FeatureMatrix(
#   fragments = f1,
#   features = gc,
#   cells = rownames(s1)
# )
# 
# p2<-FeatureMatrix(
#   fragments = f2,
#   features = gc,
#   cells = rownames(s2)
# )
# 
# 
# c1<-CreateChromatinAssay(counts = p1, sep = c(":", "-"),fragments = f1)
# counts1<-CreateSeuratObject(c1, assay = "ATAC",meta.data = s1)
# 
# 
# 
# c2<-CreateChromatinAssay(counts = p2, sep = c(":", "-"),fragments = f2)
# counts2<-CreateSeuratObject(c2, assay = "ATAC",meta.data = s2)
# 
# counts1$dataset <- 'WT'
# counts2$dataset <- 'KD'
# 
# 
# GetAssayData(counts1,slot="counts")[1:5,1:5]
# GetAssayData(counts2,slot="counts")[1:5,1:5]
# 
# combined <- merge(
#   x = counts1,
#   y = counts2,
#   add.cell.ids = c("WT","KD")
# )
# 
# combined <- RunTFIDF(combined)
# combined <- FindTopFeatures(combined, min.cutoff = 20)
# combined <- RunSVD(combined)
# combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
# 
# pdf("merge-by-pmat.Seurat.pdf")
# DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
# dev.off()
# 
# saveRDS(combined,"Merge-by-FeatureMatrix.rds")
```

![image](https://user-images.githubusercontent.com/55969398/117751292-6a410400-b247-11eb-88b5-cc61ac28af7c.png)


### Self-Merge 

```
#-----------------#-----------------
# Self-Merge:: filter both samples by gene
#   - 1. Pmat:: 1. make a csv file for scie-2-gene anno.
#   - 2. Pmat:: 2. scie-2-gene:: annotation
#   - 3. Pmat:: 3. filter pmat by annotated scie-2-g file
#   - 4. Pmat:: 4:combine/join/merge
#   - 5. singlecell.csv:: (not filtered: just rbind, will filter later in Seurat Obj.)
#   - 6. fragments.tsv (filter peaks)
#-----------------#-----------------

#-----------------#-----------------
# Pmat:: 1. make a csv file for scie-2-gene anno.
#-----------------#-----------------
# rm(list=ls())
# library(data.table)
# library(Seurat)
# library(GenomeInfoDb)
# library(EnsDb.Hsapiens.v75)
# library(ggplot2)
# library(patchwork)
# library(stringr)
# 
# x.sp1<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/L1_x.sp_pmat.snap.rds")
# 
# 
# peak1<-x.sp1@peak@elementMetadata$name
# a1<-str_replace(peak1,"b'","")
# peak1<-str_replace(a1,"'","")
# 
# chr<-mclapply(peak1,function(x){
#   str_split(x,":")[[1]][1]
# })
# chr<-unlist(chr)
# 
# len<-mclapply(peak1,function(x){
#   str_split(x,":")[[1]][2]
# })
# len<-unlist(len)
# 
# start<-mclapply(len,function(x){
#   str_split(x,"-")[[1]][1]
# })
# start<-unlist(start)
# 
# end<-mclapply(len,function(x){
#   str_split(x,"-")[[1]][2]
# })
# end<-unlist(end)
# 
# df<-data.frame(
#   chr=chr,
#   start=start,
#   end=end
# )
# rownames(df)<-peak1
# 
# write.csv(df,"L1-peak-df.csv")
# print("L1-finished")
#   
#   
# rm(list=ls())
# x.sp2<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/L2_x.sp_pmat.snap.rds")
# peak2<-x.sp2@peak@elementMetadata$name
# a2<-str_replace(peak2,"b'","")
# peak2<-str_replace(a2,"'","")
# 
# chr<-mclapply(peak2,function(x){
#   str_split(x,":")[[1]][1]
# })
# chr<-unlist(chr)
# 
# len<-mclapply(peak2,function(x){
#   str_split(x,":")[[1]][2]
# })
# len<-unlist(len)
# 
# start<-mclapply(len,function(x){
#   str_split(x,"-")[[1]][1]
# })
# start<-unlist(start)
# 
# end<-mclapply(len,function(x){
#   str_split(x,"-")[[1]][2]
# })
# end<-unlist(end)
# 
# df<-data.frame(
#   chr=chr,
#   start=start,
#   end=end
# )
# rownames(df)<-peak2
# write.csv(df2,"L2-peak-df.csv")
# print("L2 finish")


#-----------------#-----------------
# Pmat:: 2. scie-2-gene:: annotation
#-----------------#-----------------


## Python:: python3.8 peak_anno.py L1-peak-df.csv
# 
# from scie2g.csv import CSV
# 
# def peak(name):
#   bed=Csv(name,"chr","start","end")
# bed.set_annotation_from_file('hsapiens_gene_ensembl-GRCh37.p13.csv')
# bed.assign_locations_to_genes()
# bed.save_loc_to_csv(name+".csv")
# if __name__ == "__main__":
#   import sys
# peak(sys.argv[1])



#-----------------#-----------------
# Pmat:: 3. filter pmat by annotated scie-2-g file
#-----------------#-----------------

# rm(list=ls())
# library(stringr)
# 
# 
# p1<-read.csv("/Users/apple/Desktop/from-server/L1-peak-gene-anno.csv")
# #p1<-read.csv("/media/WorkingSpace/Wenyu/cardiomyocyte/5_R/Merge/CSV/L1-peak-gene-anno.csv")
# 
# p1$seq<-str_c(p1$chr,":",p1$start,"-",p1$end-1)
# 
# uniq_peaks<-c()
# for(i in unique(p1$gene_idx)){
#   peaks<-p1[p1$gene_idx==i,]
#   chr_c<-(peaks$end-peaks$start)/2
#   gene_c<-(peaks$end_position-peaks$start_position)/2
#   dis<-abs(gene_c-chr_c)
#   names(dis)<-peaks$seq
#   peaknames<-names(dis[dis==min(dis)])
#   uniq_peaks<-append(uniq_peaks,peaknames)
# }
# 
# nrow(p1)
# length(uniq_peaks)
# 
# pmat<-readRDS("L1-pmat.rds")
# pfiltered<-pmat[which(uniq_peaks %in% rownames(pmat)),]
# nrow(pfiltered)
# 
# saveRDS(pfiltered,"L1-pmat-filtered.rds")

  
# rm(list=ls())
# p2<-read.csv("/Users/apple/Desktop/from-server/L2-peak-gene-anno.csv")
# #p2<-read.csv("/media/WorkingSpace/Wenyu/cardiomyocyte/5_R/Merge/CSV/L2-peak-gene-anno.csv")
# 
# 
# p2$seq<-str_c(p2$chr,":",p2$start,"-",p2$end-1)
# 
# uniq_peaks2<-c()
# for(i in unique(p2$gene_idx)){
#   peaks<-p2[p2$gene_idx==i,]
#   chr_c<-(peaks$end-peaks$start)/2
#   gene_c<-(peaks$end_position-peaks$start_position)/2
#   dis<-abs(gene_c-chr_c)
#   names(dis)<-peaks$seq
#   peaknames<-names(dis[dis==min(dis)])
#   uniq_peaks2<-append(uniq_peaks2,peaknames)
# }
# 
# pmat2<-readRDS("L2-pmat.rds")
# pfiltered2<-pmat2[which(uniq_peaks2 %in% rownames(pmat2)),]
# nrow(pfiltered2)
# saveRDS(pfiltered2,"L2-pmat-filtered.rds")


#-----------------#-----------------
# Pmat:: 4:combine/join/merge
#-----------------#-----------------
# rm(list=ls())
# setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/pmat")
# 
# library(stringr)
# p1<-readRDS("L1-pmat-filtered.rds")
# p2<-readRDS("L2-pmat-filtered.rds")
# 
# df1<-as.data.frame(as.matrix(p1))
# colnames(df1)<-str_c("WT_",colnames(df1))
# df1$chr<-rownames(df1)
# 
# df2<-as.data.frame(as.matrix(p2))
# colnames(df2)<-str_c("KD_",colnames(df2))
# df2$chr<-rownames(df2)
# 
# 
# dmerge<-plyr::join(df1,df2,by="chr",type="full")
# rownames(dmerge)<-dmerge$chr
# dmerge$chr<-NULL
# 
# dmtx<-as.matrix(dmerge)
# dmtx<-as(dmtx,"dgCMatrix")
# 
# saveRDS(dmtx,"Merged-Pmat.rds")
# 

#-----------------#-----------------
# singlecell.csv:: (not filtered: just rbind, will filter later in Seurat Obj.)
#-----------------#-----------------
# rm(list=ls())
# library(stringr)
# setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte")
# 
# s1 <- read.csv(
#   file = "Library_20190731_001_singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )
#   
# bar1<-s1[1,]
# bar1$cell_id<-0
# 
# rownames(s1)<-str_c("WT_",rownames(s1))
# 
# 
# s2 <- read.csv(
#   file = "Library_20190731_002_singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )
# 
# bar2<-s2[1,]
# bar2$cell_id<-0
# 
# rownames(s2)<-str_c("KD_",rownames(s2))
# 
# s1<-s1[-c(1),]
# s2<-s2[-c(1),]
# 
# smerge<-rbind(s1,s2)
# 
# barcodes<-bar1+bar2
# barcodes$cell_id<-"None"
# 
# smerge<-rbind(barcodes,smerge)
# 
# 
# write.csv(smerge,"Merged-single-cell.csv")

#-----------------#-----------------
# fragments.tsv (filter peaks)
#-----------------#-----------------
# rm(list=ls())
# setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/fragments")
# library(readr)
# library(stringr)
# f1<-read_tsv("Library_20190731_001_fragments.tsv",col_names = FALSE)
# f2<-read_tsv("Library_20190731_002_fragments.tsv",col_names = FALSE)
# f1<-as.data.frame(f1)
# f2<-as.data.frame(f2)
# 
# f1$X4<-str_c("WT_",f1$X4)
# f2$X4<-str_c("KD_",f2$X4)
# 
# f1$seq<-str_c(f1$X1,":",f1$X2,"-",f1$X3)
# f2$seq<-str_c(f2$X1,":",f2$X2,"-",f2$X3)
# 
# fmerge<-rbind(f1,f2)
# 
# pmerge<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/pmat/Merged-Pmat.rds")
# peaks<-rownames(pmerge)


# seqnames<-fmerge$seq
# sum(seqnames %in% peaks)
# 
# fnew<-fmerge[which(seqnames %in% peaks),]
# write_tsv(fmerge,"Merged-fragments.tsv")
```

![Combined-dataset](https://user-images.githubusercontent.com/55969398/117751257-5b5a5180-b247-11eb-8bcb-535659da9f12.png)


### Seurat:: Merge from two filtered pmat (seems not right: ignored)

![image](https://user-images.githubusercontent.com/55969398/117751335-804ec480-b247-11eb-9e84-06eefc212ecd.png)


# Integration of Merged scATAC

![Combined-Integrated-rna+atac predicted](https://user-images.githubusercontent.com/55969398/117751454-bb50f800-b247-11eb-9bb7-fb360c680eaa.png)

![Combined-Integrated-rna+atac-by-sample](https://user-images.githubusercontent.com/55969398/117751470-c3109c80-b247-11eb-9a31-726b3843ab69.png)

![Combined-Integrated-atac predicted-by-sample](https://user-images.githubusercontent.com/55969398/117751489-cb68d780-b247-11eb-919a-423008600a89.png)

![Combined-Integrated-atac-atac predicted](https://user-images.githubusercontent.com/55969398/117751501-cdcb3180-b247-11eb-8748-1c1f01367322.png)

![Combined-TSS](https://user-images.githubusercontent.com/55969398/117751512-d3287c00-b247-11eb-9956-d378665bf1fd.png)
![Combined-Nucleosome](https://user-images.githubusercontent.com/55969398/117751517-d4f23f80-b247-11eb-86ba-4a9bd88094ef.png)
![Combined-DepthCor(atac)](https://user-images.githubusercontent.com/55969398/117751519-d6236c80-b247-11eb-980d-6d9e4f9406cd.png)
![Combined-cluster](https://user-images.githubusercontent.com/55969398/117751527-d885c680-b247-11eb-9b52-96c12ad45e69.png)
![Combined-cluster-](https://user-images.githubusercontent.com/55969398/117751533-da4f8a00-b247-11eb-9a33-ad51c618e1ff.png)
![Combined-cluster--](https://user-images.githubusercontent.com/55969398/117751538-dc194d80-b247-11eb-8e15-7298b759f9f9.png)
![Combined-c('TSS enrichment','nucleosome_signal')](https://user-images.githubusercontent.com/55969398/117751546-dde31100-b247-11eb-98f2-052cd63be8b8.png)
![Combined-c('pct_reads_in_peaks', 'peak_region_fragments')](https://user-images.githubusercontent.com/55969398/117751555-dfacd480-b247-11eb-9358-2b3cc5ac47d3.png)
![Combined-c('blacklist_ratio')](https://user-images.githubusercontent.com/55969398/117751559-e1769800-b247-11eb-9b41-a6cd0f330f34.png)


```
#-----------------#-----------------
# Integration of scATAC scRNA for merged Obj.
#   - 2. Merged by FeatureMatrix
#-----------------#-----------------
#-----------------#
# 2. Merged by FeatureMatrix
#-----------------#
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(stringr)

atac<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/Combined-atac.Seurat.rds")
DimPlot(atac,group.by = "dataset")
atac[["ATAC"]] #hromatinAssay data with 348693 features for 107234 cells

atac[["sample"]]=atac[["dataset"]]
atac[["day"]]=9
atac[["orig.ident"]]="scATAC"


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
Annotation(atac) <- annotations

atac <- NucleosomeSignal(object = atac)
atac <- TSSEnrichment(object = atac, fast = FALSE)


atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

atac$high.tss <- ifelse(atac$TSS.enrichment > 1.4, 'High', 'Low')
TSSPlot(atac, group.by = 'high.tss') + NoLegend()

atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = atac, group.by = 'nucleosome_group')

VlnPlot(
  object = atac,
  features = c('TSS.enrichment','nucleosome_signal'),
  pt.size = 0.1,
  ncol = 2,group.by = "dataset"
)

VlnPlot(
  object = atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments'),
  ncol = 2,group.by = "dataset"
)

VlnPlot(
  object = atac,
  features = c('blacklist_ratio'),
  pt.size = 0.1,
  ncol = 1,group.by = "dataset"
)

atac # 348693 features across 107234 samples within 1 assay 

atac <- subset(
  x = atac,
  subset = peak_region_fragments > 100 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.001 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1.0
)

atac #348693 features across 32919 samples within 1 assay 

atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac)

DepthCor(atac)

atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 1:25)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 1:25)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)
cluster<-Idents(atac)
table(cluster)

# 0     1     2     3     4     5     6     7     8     9    10    11 
# 12935  7787  2038  1723  1530  1448  1376  1142   983   844   678   435 

DimPlot(object = atac, label = TRUE,label.box = TRUE) 
DimPlot(object = atac, group.by = "sample")
DimPlot(object = atac, split.by = "sample") 

gene.activities <- GeneActivity(atac)
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)

DefaultAssay(atac) <- 'RNA'


rna<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/scRNA/scrna-in-one-Seurat-day-7-9.rds")
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  reduction = 'cca',
  k.anchor = 25,
  k.filter=500,
  k.score=50,
  reference.assay="RNA",
  query.assay = "RNA",
  dims=1:25
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$celltype,
  weight.reduction = atac[['lsi']],
  dims = 1:30
)


atac <- AddMetaData(object = atac, metadata = predicted.labels)
saveRDS(atac,"Combined-Integrated-by-FeatureMatrix.Seurat.rds")
```
