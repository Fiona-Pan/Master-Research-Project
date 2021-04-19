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

![barcode-selection](https://user-images.githubusercontent.com/55969398/115208418-1e5edb80-a12f-11eb-842c-d0e0c4caaa6d.png =250x)

- number of barcodes: 7010
- number of bins: 535551  #cell-by-bin matrix
- number of genes: 0
- number of peaks: 0
- number of motifs: 0

Determine significant components:

- dims = 6

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




