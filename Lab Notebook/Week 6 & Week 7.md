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

barcode selection:

```
barcodes.sel=barcodes1[which(UMI>=3.75 & UMI <=5.25 & promoter_ratio >=0.15 & promoter_ratio <=0.4),]
```

![barcode-selection](https://user-images.githubusercontent.com/55969398/115208418-1e5edb80-a12f-11eb-842c-d0e0c4caaa6d.png)



Determine significant components:

- dims = 6

![Elbow](https://user-images.githubusercontent.com/55969398/115211177-e4db9f80-a131-11eb-947b-535309b325c3.png)
![Elbow2](https://user-images.githubusercontent.com/55969398/115211676-603d5100-a132-11eb-94b7-c27f78ba036c.png)

