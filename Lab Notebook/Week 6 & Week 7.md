# Week 6 & Week 7

- apply previous scATAC-seq analysis on cardiomyocyte dataset

**Dataset**

- developmental cardiomyocyte in vitro scATAC-seq
- two sample:  
  - 1- control &emsp; 2- TMEM88 knockdown

**Methods**

- pre-processing of raw scATAC-seq data
- feature matrix construction
- dimensionality reduction
- clustering and visualization
- evaluation of reproducibility of peaks

**Results**

#### Sample1


raw barcodes:

- Total  number of barcodes: 99884

- Median number of sequencing fragments: 817

- Median number of uniquely mapped fragments: 485

- Median number of mappability ratio: 0.84

- Median number of properly paired ratio: 1

- Median number of duplicate ratio: 0.26

- Median number of chrM ratio: 0.21

- Median number of unique molecules (UMI): 485


barcode selection:

```
barcodes.sel=barcodes1[which(UMI>=3 & UMI <=5.5 & promoter_ratio >=0.1 & promoter_ratio <=0.5),]

```

![barcode-selection](https://user-images.githubusercontent.com/55969398/115208418-1e5edb80-a12f-11eb-842c-d0e0c4caaa6d.png)

- number of barcodes: 7010
- number of bins: 535551  #cell-by-bin matrix
- number of genes: 0
- number of peaks: 0
- number of motifs: 0

Determine significant components:

- dims = 8

![Elbow](https://user-images.githubusercontent.com/55969398/115211177-e4db9f80-a131-11eb-947b-535309b325c3.png)
![Dims](https://user-images.githubusercontent.com/55969398/115217323-f9bb3180-a137-11eb-867b-d63162303a01.png)

Clustering:

- UMAP:

![UMAP](https://user-images.githubusercontent.com/55969398/115217598-3e46cd00-a138-11eb-9160-3117246dd2bc.png)
![UMAP-feature](https://user-images.githubusercontent.com/55969398/115217749-633b4000-a138-11eb-9e90-51639d97b320.png)

- t-SNE:

![t-SNE](https://user-images.githubusercontent.com/55969398/115218020-b1504380-a138-11eb-8300-ef839b9567f3.png)
![t-SNE-feature](https://user-images.githubusercontent.com/55969398/115218224-e492d280-a138-11eb-8d8e-81eb9ba25dca.png)

- dendrogram:

![dendrogram](https://user-images.githubusercontent.com/55969398/115218406-1b68e880-a139-11eb-9cf3-95f3bc4fadca.png)

```
# number of cells in each cluster
  1   2   3   4   5   6   7   8   9  10  11 
807 782 739 718 687 668 641 582 541 470 375 
```



#### Sample2
raw barcodes:

- Total  number of barcodes: 109841
- Median number of sequencing fragments: 541
- Median number of uniquely mapped fragments: 269
- Median number of mappability ratio: 0.78
- Median number of properly paired ratio: 1
- Median number of duplicate ratio: 0.3
- Median number of chrM ratio: 0.44
- Median number of unique molecules (UMI): 269

barcode selection:

```
barcodes.sel=barcodes1[which(UMI>=3.5 & UMI <=5.75 & promoter_ratio >=0.1 & promoter_ratio <=0.5),]

```

![barcode](https://user-images.githubusercontent.com/55969398/115235746-1b72e380-a14d-11eb-962b-e71f3aa44796.png)

- number of barcodes: 4830
- number of bins: 535618
- number of genes: 0
- number of peaks: 0
- number of motifs: 0


Determine significant components:

- dims=10

![Elbow](https://user-images.githubusercontent.com/55969398/115235901-49582800-a14d-11eb-8eb8-35eb5f51d1d5.png)
![dims](https://user-images.githubusercontent.com/55969398/115235917-4c531880-a14d-11eb-994c-cc67fdb835a0.png)

Clustering:

- UMAP:

![UMAP](https://user-images.githubusercontent.com/55969398/115235990-668cf680-a14d-11eb-9941-68ed4ad600cd.png)
![UMAP-feature](https://user-images.githubusercontent.com/55969398/115236028-73114f00-a14d-11eb-85e6-a63d79466229.png)


- t-SNE:

![t-SNE](https://user-images.githubusercontent.com/55969398/115236083-86241f00-a14d-11eb-9817-2cf406e18f8e.png)
![t-SNE-feature](https://user-images.githubusercontent.com/55969398/115236139-94723b00-a14d-11eb-9372-57760cbeb2bf.png)

- dendrogram:

![dendrogram](https://user-images.githubusercontent.com/55969398/115253268-cb9d1800-a15e-11eb-9a2e-9563417d5994.png)

```
  1   2   3   4   5   6   7   8   9  10  11  12  13 
608 571 513 476 426 400 395 358 250 235 226 191 181 
```

cell-by-peak matrix:

- dimension: 4830(cells) x 238842 (peaks)
- non-zero count: 0.07618712
- by-cluster:

```
[1] "non-zero peaks couts: cluster1: 0.0588"
[1] "non-zero peaks couts: cluster2: 0.067"
[1] "non-zero peaks couts: cluster3: 0.1588"
[1] "non-zero peaks couts: cluster4: 0.1035"
[1] "non-zero peaks couts: cluster5: 0.0803"
[1] "non-zero peaks couts: cluster6: 0.0633"
[1] "non-zero peaks couts: cluster7: 0.069"
[1] "non-zero peaks couts: cluster8: 0.0623"
[1] "non-zero peaks couts: cluster9: 0.0734"
[1] "non-zero peaks couts: cluster10: 0.0199"
[1] "non-zero peaks couts: cluster11: 0.0534"
[1] "non-zero peaks couts: cluster12: 0.0306"
[1] "non-zero peaks couts: cluster13: 0.0731"
```


#### Selection of high-value peaks in each cluster

```
#-----------------
# Extract high-score peaks
#-----------------
library(ChIPseeker)

setwd("/Users/apple/Downloads/Master_Research/data/Cardiomyocyte/L1peak/edited")
peaks.names = system("ls *narrowPeak| grep narrowPeak", intern=TRUE)
names=str_match(peaks.names,"Library_20190731_001_S2.(.*)_peaks.edited.narrowPeak")[,2]
peak.gr.ls = lapply(peaks.names, function(x){ # a list of GRanges of clusters
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]),names=peak.df[,4],score=peak.df[,5],signalValue=peak.df[,7],pvalue=peak.df[,8],qvalue=peak.df[,9],peak=peak.df[,10])
})
names(peak.gr.ls)<-names
# how many peaks for each cluster
for(i in 1:length(peak.gr.ls)){
  print(length(peak.gr.ls[[i]]))
}

# quantile for each cluster - by score
for(i in 1:length(peak.gr.ls)){
  print(quantile(peak.gr.ls[[i]]$score))
} # selection of 75% as threshold


# quantile for each cluster - by signalVal
for(i in 1:length(peak.gr.ls)){
  print(quantile(peak.gr.ls[[i]]$signalValue))
} # selection of 50% as threshold


# new list of high-score cluster peaks
q=list()
for(i in 1:length(peak.gr.ls)){
  new=peak.gr.ls[[i]][peak.gr.ls[[i]]$score>quantile(peak.gr.ls[[i]]$score)["75%"]]
  q=append(q,list(new))
}

# new list of high-signalVal cluster peaks
s=list()
for(i in 1:length(q)){
  new=q[[i]][q[[i]]$signalValue>quantile(q[[i]]$signalValue)["50%"]]
  s=append(s,list(new))
}

# check
for(i in 1:length(s)){
  print(length(s[[i]]))
}

# peak.gr.ls.new
peak.gr.ls.new=s
names(peak.gr.ls.new)=names

# peaks combined
peak.gr = reduce(Reduce(c, peak.gr.ls.new)) # reduce into one GRange object for all clusters
peak.gr

# save each cluster peaks as bed
for(i in 1:length(peak.gr.ls.new)){
  peak.gr.df=as.data.frame(peak.gr.ls.new[[i]])
  peak.gr.df[,"strand"]='.'
  peak.gr.df=peak.gr.df[,c(1,2,3,6,7,5,8,9,10,11)]
  write.table(peak.gr.df,file=str_c("peaks",names[i],".narrowPeak"),append=FALSE,
              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".",
              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
}

#save combined bed
peaks.df = as.data.frame(peak.gr)[,1:3]
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".",
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

# snaptools snap-add-pmat \
# --snap-file atac_v1_pbmc_5k.snap \
# --peak-file peaks.combined.bed 


# peaks plot
peak10=peak.gr.ls.new[[1]]
covplot(peak10,weightCol="signalValue")
covplot(peak.gr.ls.new,weightCol="signalValue")

```


