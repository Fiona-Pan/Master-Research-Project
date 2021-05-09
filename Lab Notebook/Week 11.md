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
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
#-----------------
plan("multiprocess", workers = 2)

p1<-read.csv("/Users/apple/Desktop/from-server/CSV/L1-peak-df.csv")
p1$X<-NULL
p1$C1<-NULL
p1$C2<-NULL

g1<-makeGRangesFromDataFrame(p1)
#-----------------
p2<-read.csv("/Users/apple/Desktop/from-server/CSV/L2-peak-df.csv")
p2$X<-NULL
p2$C1<-NULL
p2$C2<-NULL

g2<-makeGRangesFromDataFrame(p2)
#-----------------
combined.peaks <- reduce(x = c(g1,g2))
#-----------------
s1<-read.table(
  file = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_001_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
dim(s1)
s1<-s1[s1$passed_filters > 100, ]

#-----------------
s2<-read.table(
  file = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_002_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]
dim(s2)
s2<-s2[s2$passed_filters > 100, ]

#-----------------
f1 <- CreateFragmentObject(
  path = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_001_fragments.tsv.gz",
  cells = rownames(s1)
)
#-----------------
f2 <- CreateFragmentObject(
  path = "/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/Library_20190731_002_fragments.tsv.gz",
  cells = rownames(s2)
)
#-----------------
counts1<-FeatureMatrix(
  fragments = f1,
  features = combined.peaks,
  cells = rownames(s1)
)

dim(counts1) # 348693  68435
#-----------------
counts2<-FeatureMatrix(
  fragments = f2,
  features = combined.peaks,
  cells = rownames(s2)
)

dim(counts2)


#-----------------
atac1_assay <- CreateChromatinAssay(counts1, fragments = f1)
atac1 <- CreateSeuratObject(atac1_assay, assay = "ATAC")

# 348693 features across 68435 samples within 1 assay 


#-----------------
atac2_assay <- CreateChromatinAssay(counts2, fragments = f2)
atac2 <- CreateSeuratObject(atac2_assay, assay = "ATAC")

# 348693 features across 38799 samples within 1 assay 

#-----------------

atac1$dataset <- 'WT'
atac2$dataset <- 'KD'

combined <- merge(
  x = atac1,
  y = atac2,
  add.cell.ids = c("WT","KD")
)
# 348693 features across 107234 samples within 1 assay 

combined[["ATAC"]]
#-----------------
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')


DimPlot(combined, group.by = 'dataset')

CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr1-569860-570219"
)
#-----------------
saveRDS(combined,"Combined-atac.Seurat.rds")
```


### Merge by-hand

```
# -----------------#-----------------
# Pmat:: 1. make a csv file for scie-2-gene anno.
# -----------------#-----------------
rm(list=ls())
library(data.table)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(stringr)

x.sp1<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/L1_x.sp_pmat.snap.rds")


peak1<-x.sp1@peak@elementMetadata$name
a1<-str_replace(peak1,"b'","")
peak1<-str_replace(a1,"'","")

chr<-mclapply(peak1,function(x){
  str_split(x,":")[[1]][1]
})
chr<-unlist(chr)

len<-mclapply(peak1,function(x){
  str_split(x,":")[[1]][2]
})
len<-unlist(len)

start<-mclapply(len,function(x){
  str_split(x,"-")[[1]][1]
})
start<-unlist(start)

end<-mclapply(len,function(x){
  str_split(x,"-")[[1]][2]
})
end<-unlist(end)

df<-data.frame(
  chr=chr,
  start=start,
  end=end
)
rownames(df)<-peak1

write.csv(df,"L1-peak-df.csv")
print("L1-finished")
#-----------------#-----------------
#-----------------#-----------------
rm(list=ls())
x.sp2<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/L2_x.sp_pmat.snap.rds")
peak2<-x.sp2@peak@elementMetadata$name
a2<-str_replace(peak2,"b'","")
peak2<-str_replace(a2,"'","")

chr<-mclapply(peak2,function(x){
  str_split(x,":")[[1]][1]
})
chr<-unlist(chr)

len<-mclapply(peak2,function(x){
  str_split(x,":")[[1]][2]
})
len<-unlist(len)

start<-mclapply(len,function(x){
  str_split(x,"-")[[1]][1]
})
start<-unlist(start)

end<-mclapply(len,function(x){
  str_split(x,"-")[[1]][2]
})
end<-unlist(end)

df<-data.frame(
  chr=chr,
  start=start,
  end=end
)
rownames(df)<-peak2
write.csv(df2,"L2-peak-df.csv")
print("L2 finish")
# -----------------#-----------------
# Pmat:: 2. scie-2-gene:: annotation
# -----------------#-----------------
# Python:: python3.8 peak_anno.py L1-peak-df.csv

from scie2g.csv import CSV

def peak(name):
  bed=Csv(name,"chr","start","end")
bed.set_annotation_from_file('hsapiens_gene_ensembl-GRCh37.p13.csv')
bed.assign_locations_to_genes()
bed.save_loc_to_csv(name+".csv")
if __name__ == "__main__":
  import sys
peak(sys.argv[1])
# -----------------#-----------------
# Pmat:: 3. filter pmat by annotated scie-2-g file
# -----------------#-----------------
rm(list=ls())
library(stringr)


p1<-read.csv("/Users/apple/Desktop/from-server/L1-peak-gene-anno.csv")
#p1<-read.csv("/media/WorkingSpace/Wenyu/cardiomyocyte/5_R/Merge/CSV/L1-peak-gene-anno.csv")

p1$seq<-str_c(p1$chr,":",p1$start,"-",p1$end-1)

uniq_peaks<-c()
for(i in unique(p1$gene_idx)){
  peaks<-p1[p1$gene_idx==i,]
  chr_c<-(peaks$end-peaks$start)/2
  gene_c<-(peaks$end_position-peaks$start_position)/2
  dis<-abs(gene_c-chr_c)
  names(dis)<-peaks$seq
  peaknames<-names(dis[dis==min(dis)])
  uniq_peaks<-append(uniq_peaks,peaknames)
}

nrow(p1)
length(uniq_peaks)

pmat<-readRDS("L1-pmat.rds")
pfiltered<-pmat[which(uniq_peaks %in% rownames(pmat)),]
nrow(pfiltered)

saveRDS(pfiltered,"L1-pmat-filtered.rds")
# -----------------#-----------------
rm(list=ls())
p2<-read.csv("/Users/apple/Desktop/from-server/L2-peak-gene-anno.csv")
#p2<-read.csv("/media/WorkingSpace/Wenyu/cardiomyocyte/5_R/Merge/CSV/L2-peak-gene-anno.csv")


p2$seq<-str_c(p2$chr,":",p2$start,"-",p2$end-1)

uniq_peaks2<-c()
for(i in unique(p2$gene_idx)){
  peaks<-p2[p2$gene_idx==i,]
  chr_c<-(peaks$end-peaks$start)/2
  gene_c<-(peaks$end_position-peaks$start_position)/2
  dis<-abs(gene_c-chr_c)
  names(dis)<-peaks$seq
  peaknames<-names(dis[dis==min(dis)])
  uniq_peaks2<-append(uniq_peaks2,peaknames)
}

pmat2<-readRDS("L2-pmat.rds")
pfiltered2<-pmat2[which(uniq_peaks2 %in% rownames(pmat2)),]
nrow(pfiltered2)
saveRDS(pfiltered2,"L2-pmat-filtered.rds")
# -----------------#-----------------
# Pmat:: 4:combine/join/merge
# -----------------#-----------------
rm(list=ls())
setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/pmat")

library(stringr)
p1<-readRDS("L1-pmat-filtered.rds")
p2<-readRDS("L2-pmat-filtered.rds")

df1<-as.data.frame(as.matrix(p1))
colnames(df1)<-str_c("WT_",colnames(df1))
df1$chr<-rownames(df1)

df2<-as.data.frame(as.matrix(p2))
colnames(df2)<-str_c("KD_",colnames(df2))
df2$chr<-rownames(df2)


dmerge<-plyr::join(df1,df2,by="chr",type="full")
rownames(dmerge)<-dmerge$chr
dmerge$chr<-NULL

dmtx<-as.matrix(dmerge)
dmtx<-as(dmtx,"dgCMatrix")

saveRDS(dmtx,"Merged-Pmat.rds")

# -----------------
# singlecell.csv:: (not filtered: just rbind, will filter later in Seurat Obj.)
# -----------------
rm(list=ls())
library(stringr)
setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte")

s1 <- read.csv(
  file = "Library_20190731_001_singlecell.csv",
  header = TRUE,
  row.names = 1
)

bar1<-s1[1,]
bar1$cell_id<-0

rownames(s1)<-str_c("WT_",rownames(s1))


s2 <- read.csv(
  file = "Library_20190731_002_singlecell.csv",
  header = TRUE,
  row.names = 1
)

bar2<-s2[1,]
bar2$cell_id<-0

rownames(s2)<-str_c("KD_",rownames(s2))

s1<-s1[-c(1),]
s2<-s2[-c(1),]

smerge<-rbind(s1,s2)

barcodes<-bar1+bar2
barcodes$cell_id<-"None"

smerge<-rbind(barcodes,smerge)


write.csv(smerge,"Merged-single-cell.csv")
# # -# -# -# -# -# -# -# -# -# -# -# -# -# -# -# -# -
# fragments.tsv (filter peaks)
# # -# -# -# -# -# -# -# -# -# -# -# -# -# -# -# -# -
rm(list=ls())
setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/fragments")
library(readr)
library(stringr)
f1<-read_tsv("Library_20190731_001_fragments.tsv",col_names = FALSE)
f2<-read_tsv("Library_20190731_002_fragments.tsv",col_names = FALSE)
f1<-as.data.frame(f1)
f2<-as.data.frame(f2)

f1$X4<-str_c("WT_",f1$X4)
f2$X4<-str_c("KD_",f2$X4)

f1$seq<-str_c(f1$X1,":",f1$X2,"-",f1$X3)
f2$seq<-str_c(f2$X1,":",f2$X2,"-",f2$X3)

fmerge<-rbind(f1,f2)

pmerge<-readRDS("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/RDS/pmat/Merged-Pmat.rds")
peaks<-rownames(pmerge)
seqnames<-fmerge$seq
sum(seqnames %in% peaks)

fnew<-fmerge[which(seqnames %in% peaks),]
write_tsv(fmerge,"Merged-fragments.tsv")
```
