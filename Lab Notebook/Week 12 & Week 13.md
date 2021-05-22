# WEEK 12 & WEEK 13

### Deep Analysis based on Merged, filtered, Integrated scATAC-seq dataset

- only select peaks at promoter region (~15% genomic regions)
- DAR
- GO
- KEGG
- motif
- Main:
  - use of gene-by-cell matrix for Merged and integrated scATAC-seq dataset, split by samples (WT & KD)
  - extract genes that are differentially expressed in Cardiac Cells (vs. the rest) for two samples
  - selection of biological important (& less by chance) motifs by filtering motifs based on prediction of motifs from scRNA gene-by-cell matrix
  - use those gene lists as inputs for GO and KEGG

### selection of peaks at promoter regions only (& Cardiac cells only)

![image](https://user-images.githubusercontent.com/55969398/119228195-792b8e80-bb44-11eb-98a1-a172c523ccc5.png)

```
#-----------------#
# 1. select promoter/enhancer only (cardiac cells only)
#-----------------#
library(ChIPpeakAnno)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
atac<-readRDS("/Users/apple/Desktop/New/anchor_10_filter_NA_score_30-Final.rds")
for(i in levels(atac)) {
  cells_to_reid <- WhichCells(atac, idents = i)
  newid <- names(sort(table(atac$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(atac, cells = cells_to_reid) <- newid
}
atac$seurat_annotations<-Idents(atac)
cardiac<-subset(atac,subset=seurat_annotations=="Cardiac cells")
DefaultAssay(cardiac)<-"ATAC"
aCR<-assignChromosomeRegion(granges(cardiac), nucleotideLevel=FALSE, 
                            precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                            TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
aCR$percentage

annotationData <- promoters(TSS.human.GRCh37, upstream=1000, downstream=1000)
annotatedPeak <- annotatePeakInBatch(granges(cardiac), 
                                     AnnotationData=annotationData,
                                     output="overlapping")
cardiac.gr<-annotatedPeak[which(!is.na(annotatedPeak$fromOverlappingOrNearest)),]

cardiac.gr<-data.frame(cardiac.gr)

cardiac.peak<-str_c(cardiac.gr$seqnames,"-",cardiac.gr$start,"-",cardiac.gr$end)

cardiac.promoter<-cardiac[cardiac.peak,]
cardiac.promoter
DefaultAssay(cardiac.promoter)<-"ATAC"
cardiac.promoter<-SetIdent(cardiac.promoter,value="dataset")

```

### DAR & motif analysis

- 329 enriched motifs (p.val<0.05)
- Ex.plot: top 12 enriched motifs based on fold.enrichment

![image](https://user-images.githubusercontent.com/55969398/119228366-48982480-bb45-11eb-8c8e-7ac21b003cd8.png)


```
#-----------------#
# 2. DAR peak 
#-----------------#
cardiac.promoter<-ScaleData(cardiac.promoter)
marker.genes<-FindAllMarkers(cardiac.promoter,
                             slot="scale.data",
                             group.by = "dataset",
                             logfc.threshold = 0.20,
                             min.pct=0.1,
                             min.cells.feature=3,
                             test.use="bimod",
                             only.pos=TRUE)

# marker.genes.WT<-marker.genes %>% subset(cluster=="WT")
# marker.genes.KD<-marker.genes %>% subset(cluster=="KD")

# 
# par(mfrow=c(2,1),mar=c(5,12,3,3))
# top30<-marker.genes.WT %>% top_n(30,avg_diff)
# barplot(sort(setNames(top30$avg_diff,rownames(top30)),F),horiz=T,las=1,main="control")
# top30.KD<-marker.genes.KD %>% top_n(30,avg_diff)
# barplot(sort(setNames(top30.KD$avg_diff,rownames(top30.KD)),F),horiz=T,las=1,main="TMEM88 KD")

top.da.peak <- rownames(marker.genes[marker.genes$p_val < 0.01, ])


cardiac.promoter<-readRDS("/Users/apple/Desktop/New/cardiac.promoter.rds")
cardiac.promoter

#-----------------#
# 3. Motif
#-----------------#
library(JASPAR2020)
library(TFBSTools)
DefaultAssay(cardiac.promoter)<-"ATAC"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# cardiac.promoter<-cardiac.promoter[AccessiblePeaks(cardiac.promoter),]


cardiac.motif <- AddMotifs(
  object = cardiac.promoter,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pfm
)


enriched.motifs <- FindMotifs(
  object = cardiac.motif,
  features = top.da.peak
)

enriched.motifs<- enriched.motifs[enriched.motifs$pvalue<0.05,]
topenriched<-enriched.motifs %>% top_n(6,fold.enrichment)

MotifPlot(
  object = cardiac.motif,assay = "ATAC",
  motifs = rownames(topenriched)
)


cardiac.motif <- RunChromVAR(
  object = cardiac.motif,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

DefaultAssay(cardiac.motif) <- 'chromvar'

# cardiac.motif.filter<-cardiac.motif[rownames(enriched.motifs),]



differential.activity <- FindMarkers(
  object = cardiac.motif,
  ident.1 = 'WT',
  group.by="dataset",
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

differential.activity<-differential.activity[differential.activity$p_val<0.05,]

MotifPlot(
  object = cardiac.motif,
  motifs = rownames(differential.activity),
  assay = 'ATAC'
)

differential.activity.KD <- FindMarkers(
  object = cardiac.motif,
  ident.1 = 'KD',
  group.by="dataset",
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

differential.activity.KD<-differential.activity.KD[differential.activity.KD$p_val<0.05,]

```

### motif filteration:: scRNA-seq TF prediction

- 92 matched motifs b/t scATAC-seq motif enrichment analysis & scRNA-seq TF prediction
- targeting 2576 genes in scRNA-seq dataset (2925 genes in scATAC-seq gene-by-cell matrix)

```
# matched motifs
 [1] "ASCL1"   "ATF2"    "BACH1"   "BHLHE40" "CREB1"   "CREB3"   "CREB3L1" "CTCFL"   "E2F6"   
[10] "E2F7"    "EBF1"    "ELF1"    "ELK4"    "EOMES"   "ERG"     "ESR2"    "ESRRA"   "ETS1"   
[19] "ETS2"    "ETV1"    "ETV4"    "FLI1"    "FOSL2"   "GABPA"   "GLI2"    "HMBOX1"  "HNF1A"  
[28] "HOXA9"   "IKZF1"   "IRF1"    "IRF2"    "IRF3"    "IRF4"    "IRF9"    "JUN"     "KLF3"   
[37] "KLF5"    "KLF6"    "MAFK"    "MAX"     "MAZ"     "MEF2B"   "MEIS1"   "MEIS2"   "MITF"   
[46] "MNT"     "MYBL2"   "MYC"     "MYCN"    "NFKB1"   "NFKB2"   "NFYA"    "NFYB"    "NR2C2"  
[55] "NR2F1"   "NR5A1"   "NRF1"    "PAX5"    "PAX6"    "PBX3"    "PKNOX1"  "PRDM1"   "PROX1"  
[64] "RARA"    "RBPJ"    "RELA"    "RELB"    "REST"    "RUNX2"   "SNAI2"   "SP1"     "SP2"    
[73] "SP3"     "SPI1"    "SPIB"    "TBX21"   "TCF3"    "TCF4"    "TFAP2A"  "TFAP2C"  "TFDP1"  
[82] "THAP1"   "THAP11"  "TP53"    "TP73"    "USF1"    "VDR"     "ZBTB7A"  "ZEB1"    "ZNF143" 
[91] "ZNF263"  "ZNF740" 
```

```
#-----------------#
# 4. RNA-prediction (scRNA)
#-----------------#

library(dorothea)
library(viper)

rna<-readRDS("/Users/apple/Downloads/Master_Research/from-server/scrna-in-one-Seurat-day-8-9.rds")
rna
rna.cardiac<-subset(rna,subset=seurat_clusters=="Cardiac cells")
rna.cardiac
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))
rna.cardiac <- run_viper(rna.cardiac, regulon,
                         options = list(method = "scale", minsize = 4, 
                                        eset.filter = FALSE, cores = 1, 
                                        verbose = FALSE))


DefaultAssay(object = rna.cardiac) <- "dorothea"
rna.cardiac <- ScaleData(rna.cardiac)
viper_scores_df <- GetAssayData(rna.cardiac, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

Idents(rna.cardiac)<-rna.cardiac$dataset
CellsClusters <- data.frame(cell = names(Idents(rna.cardiac)), 
                            cell_type = as.character(Idents(rna.cardiac)),
                            check.names = F)

viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
highly_variable_tfs <- summarized_viper_scores %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(200, avg) %>%
  distinct(tf)
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

summarized_viper_scores_df$cell_type<-NULL
palette_length = 1000
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))
library(pheatmap)
viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 




sum(highly_variable_tfs$tf %in% enriched.motifs$motif.name)
motif.match<-highly_variable_tfs$tf[highly_variable_tfs$tf %in% enriched.motifs$motif.name]

genes.match<-dorothea_regulon_human[which(dorothea_regulon_human$tf %in% motif.match),]$target
genes.match<-unique(genes.match) #3109

sum(unique(genes.match) %in% AccessiblePeaks(rna.cardiac)) #2369
genes.target<-genes.match[unique(genes.match) %in% rownames(rna.cardiac)]


a1<-dorothea_regulon_human[which(dorothea_regulon_human$tf %in% motif.match),]
a1<-a1[which(a1$target %in% genes.target),]
a1<-data.frame(a1)

high.confidence<-a1[a1$confidence=="A",]
unique(high.confidence$tf)
unique(high.confidence$target)

```


```
# data.frame:: high-confidence predicted & matched motifs
     tf confidence    target mor
45   ATF2          A     APOC3   1
46   ATF2          A      ATF3   1
47   ATF2          A     CCND1   1
48   ATF2          A    CDKN1A   1
49   ATF2          A     DUSP1   1
50   ATF2          A       FN1   1
51   ATF2          A       FOS   1
52   ATF2          A    IMPDH2   1
53   ATF2          A       JUN   1
54   ATF2          A      MMP2   1
55   ATF2          A      NOS3  -1
56   ATF2          A     NR4A2   1
57   ATF2          A      PCK2   1
58   ATF2          A      PLAT   1
59   ATF2          A       RB1   1
60   ATF2          A     TGFB2  -1
146 CREB1          A      AQP3   1
147 CREB1          A      BCL2   1
148 CREB1          A    BDKRB1   1
149 CREB1          A      BDNF   1
150 CREB1          A     CCND1   1
151 CREB1          A       CD4   1
152 CREB1          A    CDK11A   1
153 CREB1          A     CPT1A   1
154 CREB1          A     CXCR4   1
155 CREB1          A   CYP11A1   1
156 CREB1          A      DIO2   1
157 CREB1          A      EGR1   1
158 CREB1          A     ERBB2   1
159 CREB1          A      ETV3   1
160 CREB1          A      FLT1   1
161 CREB1          A       FN1   1
162 CREB1          A       FOS   1
163 CREB1          A     HMOX1  -1
164 CREB1          A      HPGD   1
165 CREB1          A       IL6   1
166 CREB1          A       JUN  -1
167 CREB1          A      MITF   1
168 CREB1          A      MSMB   1
169 CREB1          A      MUC4   1
170 CREB1          A       NF1   1
171 CREB1          A     NOLC1   1
172 CREB1          A       NPY   1
173 CREB1          A     NR4A2   1
174 CREB1          A     NR4A3   1
175 CREB1          A       NTS   1
176 CREB1          A      ODC1   1
177 CREB1          A      PCK2   1
178 CREB1          A      PLAT  -1
179 CREB1          A     POLD2   1
180 CREB1          A    PPP2CA   1
181 CREB1          A     PSEN1   1
182 CREB1          A      RARB   1
183 CREB1          A     RRM2B   1
184 CREB1          A   SLC19A1   1
185 CREB1          A   SLC20A1   1
186 CREB1          A     SPRY2   1
187 CREB1          A      STAR   1
188 CREB1          A       TRH   1
189 CREB1          A      UGT8   1
190 CREB1          A       UXT   1
191 CREB1          A     VEGFA   1
531   ERG          A      CDH5   1
532   ERG          A       ENG   1
533   ERG          A      FGF2  -1
534   ERG          A     HMOX1   1
535   ERG          A     ICAM1  -1
536   ERG          A     ICAM2   1
537   ERG          A      MMP9   1
538   ERG          A      NOS3   1
539   ERG          A   TMPRSS2   1
540   ERG          A       VIM   1
541   ERG          A       VWF   1
542  ESR2          A      CD68   1
543  ESR2          A   COLEC12   1
544  ESR2          A      CTSD   1
545  ESR2          A     EBAG9   1
546  ESR2          A       FOS   1
547  ESR2          A    IL1RAP   1
548  ESR2          A       JUN   1
549  ESR2          A     MKNK2   1
550  ESR2          A       OXT   1
551  ESR2          A      RARA   1
552  ESR2          A  SERPINE1  -1
553  ESR2          A      TERT   1
554  ESR2          A      TFF1   1
555  ESR2          A      TGFA   1
556  ESR2          A    TMSB4X   1
618  ETS1          A     ABCC1   1
619  ETS1          A    ANGPT2   1
620  ETS1          A     ANPEP   1
621  ETS1          A     ASAP1   1
622  ETS1          A    ATP2A3   1
623  ETS1          A     ATXN2   1
624  ETS1          A    B3GAT3   1
625  ETS1          A       BAX  -1
626  ETS1          A      BCL2   1
627  ETS1          A       BID   1
628  ETS1          A     BIRC2   1
629  ETS1          A      BMP4   1
630  ETS1          A     CASP8   1
631  ETS1          A     CCND1   1
632  ETS1          A       CD4   1
633  ETS1          A     CDC37   1
634  ETS1          A     CDH13   1
635  ETS1          A    CDK11A   1
636  ETS1          A      CDK4   1
637  ETS1          A    CDKN1A   1
638  ETS1          A    CDKN1B   1
639  ETS1          A      CHUK   1
640  ETS1          A    COL1A1   1
641  ETS1          A     COPS2   1
642  ETS1          A    COX4I1   1
643  ETS1          A     CSF3R   1
644  ETS1          A   CSNK2A1   1
645  ETS1          A   CSNK2A2   1
646  ETS1          A    CSNK2B   1
647  ETS1          A      CTSB   1
648  ETS1          A     CXCR4   1
649  ETS1          A      DAD1   1
650  ETS1          A     DUSP6   1
651  ETS1          A      ECE1   1
652  ETS1          A      EGR1   1
653  ETS1          A      EGR2   1
654  ETS1          A      EMSY   1
655  ETS1          A       ERG   1
656  ETS1          A      ETS2   1
657  ETS1          A      ETV4   1
658  ETS1          A       FAS   1
659  ETS1          A      FLI1   1
660  ETS1          A      FLT1   1
661  ETS1          A       FOS   1
662  ETS1          A     FOSL1   1
663  ETS1          A     FOXD1   1
664  ETS1          A     GABPA   1
665  ETS1          A     HMOX1   1
666  ETS1          A      HPGD   1
667  ETS1          A      IER5   1
668  ETS1          A     ITGA5   1
669  ETS1          A      JAK3   1
670  ETS1          A       LCK   1
671  ETS1          A      LCP1   1
672  ETS1          A      MCL1   1
673  ETS1          A      MDM2  -1
674  ETS1          A       MET   1
675  ETS1          A     MGAT5   1
676  ETS1          A    MGAT5B   1
677  ETS1          A      MMP2   1
678  ETS1          A      MMP9   1
679  ETS1          A       MYB   1
680  ETS1          A     NDRG1  -1
681  ETS1          A     NFKB1   1
682  ETS1          A      NOS3   1
683  ETS1          A      NRP1   1
684  ETS1          A     PARP1   1
685  ETS1          A     PCSK6   1
686  ETS1          A       PF4   1
687  ETS1          A      PLAU   1
688  ETS1          A     POLD1   1
689  ETS1          A       POR   1
690  ETS1          A     PSEN1   1
691  ETS1          A     RAD51   1
692  ETS1          A     RUNX1   1
693  ETS1          A   SLC19A1   1
694  ETS1          A     SMAD4   1
695  ETS1          A     SOCS1   1
696  ETS1          A     SPRY2   1
697  ETS1          A     SURF1   1
698  ETS1          A     TAF12   1
699  ETS1          A    TBXAS1  -1
700  ETS1          A    TCEAL1   1
701  ETS1          A      TERT  -1
702  ETS1          A    TFAP2A   1
703  ETS1          A      TFRC   1
704  ETS1          A    TGFBR2   1
705  ETS1          A     THBS1   1
706  ETS1          A     TIMP1  -1
707  ETS1          A       TNC   1
708  ETS1          A TNFRSF10A   1
709  ETS1          A TNFRSF10B   1
710  ETS1          A      TOM1   1
711  ETS1          A      TP53   1
712  ETS1          A     TTYH3   1
713  ETS1          A     VEGFA   1
714  ETS1          A       VWF   1
715  ETS1          A    ZNF175   1
716  ETS2          A    ANGPT2   1
717  ETS2          A     ANPEP   1
718  ETS2          A     BRCA1   1
719  ETS2          A      CDK1   1
720  ETS2          A     CHRNE   1
721  ETS2          A     CSF3R   1
722  ETS2          A      EGR1   1
723  ETS2          A       ERG   1
724  ETS2          A      FLI1   1
725  ETS2          A       FOS   1
726  ETS2          A     FOSL1   1
727  ETS2          A     ICAM1   1
728  ETS2          A      JAK3   1
729  ETS2          A      MDM2   1
730  ETS2          A      MMP9   1
731  ETS2          A      MSR1   1
732  ETS2          A       MYC   1
733  ETS2          A     PDE7A   1
734  ETS2          A     PSEN1   1
735  ETS2          A      TERT   1
736  ETS2          A      TP53   1
737  ETS2          A       VWF   1
766  ETV4          A    BDKRB1   1
767  ETV4          A     ERBB2  -1
768  ETV4          A      ETS1   1
769  ETV4          A      MMP2   1
770  ETV4          A      MUC4   1
771  ETV4          A     PLAUR   1
772  ETV4          A      SYN2   1
773  ETV4          A    TGFBR2   1
774  ETV4          A     TIMP2   1
775  ETV4          A       VIM   1
776  FLI1          A       ERG   1
777  FLI1          A      EYA3   1
778  FLI1          A       FOS   1
779  FLI1          A     HMOX1   1
780  FLI1          A       ID2   1
781  FLI1          A    IGFBP3  -1
782  FLI1          A    ITGA2B   1
783  FLI1          A      TERT   1
784  FLI1          A    TGFBR2   1
785 FOSL2          A      BCL6   1
786 FOSL2          A     BRCA1   1
787 FOSL2          A       CLU   1
788 FOSL2          A     FOSL1   1
789 FOSL2          A      NOS2   1
790 FOSL2          A     PLAUR   1
842  GLI2          A      BCL2   1
843  GLI2          A     CCND1   1
844  GLI2          A    COL5A2   1
845  GLI2          A    EFEMP1  -1
846  GLI2          A       FAS  -1
847  GLI2          A     FGF13   1
848  GLI2          A     FOXA2   1
849  GLI2          A       FST   1
850  GLI2          A      GLI1   1
851  GLI2          A    IFNGR1  -1
852  GLI2          A       LUM  -1
853  GLI2          A      PCNA   1
854  GLI2          A     PTCH1   1
855  GLI2          A      UGCG  -1
856  GLI2          A     WNT2B   1
882 HNF1A          A       ALB   1
883 HNF1A          A     HNF4A   1
```


### GO & KEGG

![image](https://user-images.githubusercontent.com/55969398/119228792-3f0fbc00-bb47-11eb-99f3-fe43e2c8387f.png)

![image](https://user-images.githubusercontent.com/55969398/119228815-5f3f7b00-bb47-11eb-87d6-94eee203b24f.png)


```
#-----------------#
# 5. GO & KEGG
#-----------------#
DefaultAssay(cardiac.promoter)<-"RNA"
cardiac.promoter.motif.gene<-cardiac.promoter[genes.target,]
cardiac.promoter.motif.gene


marker.genes.ALL<-FindAllMarkers(cardiac.promoter.motif.gene,
                                 slot="scale.data",
                                 group.by = "dataset",
                                 logfc.threshold = 0.20,
                                 min.pct=0.1,
                                 min.cells.feature=3,
                                 test.use="bimod",
                                 only.pos=TRUE)
table(marker.genes.ALL$cluster)

marker.genes.WT<-marker.genes.ALL[which(marker.genes.ALL$cluster=="WT"),]
marker.genes.KD<-marker.genes.ALL[which(marker.genes.ALL$cluster=="KD"),]



par(mfrow=c(2,1),mar=c(5,8,3,3))
top30<-marker.genes.WT %>% top_n(30,avg_diff)
barplot(sort(setNames(top30$avg_diff,rownames(top30)),F),horiz=T,las=1,main="control")
top30.KD<-marker.genes.KD %>% top_n(30,avg_diff)
barplot(sort(setNames(top30.KD$avg_diff,rownames(top30.KD)),F),horiz=T,las=1,main="TMEM88 KD")


DoHeatmap(subset(cardiac.promoter.motif.gene,downsample=20),features=rownames(top30),size=3,angle=0,group.by = "seurat_annotations",group.colors = c("dark blue","dark blue"))+scale_fill_gradient(low="dark blue",high="light blue")+ggtitle("control")
DoHeatmap(subset(cardiac.promoter.motif.gene,downsample=30),features=rownames(top30.KD),size=3,angle=0,group.by = "seurat_annotations",group.colors = c("dark blue","dark blue"))+scale_fill_gradient(low="dark blue",high="light blue")+ggtitle("KD")


genes.WT<-rownames(marker.genes.WT)
genes.KD<-rownames(marker.genes.KD)

gene.WT<-bitr(genes.WT,fromType = "SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
gene.WT <- pull(gene.WT,ENTREZID)   

gene.KD<-bitr(genes.KD,fromType = "SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
gene.KD <- pull(gene.KD,ENTREZID) 

genelist<-list()
genelist$WT<-gene.WT
genelist$KD<-gene.KD


# KEGG
c.KEGG <- compareCluster(genelist, fun="enrichKEGG",organism="hsa", pvalueCutoff=0.05)
dotplot(c.KEGG, showCategory=20)

# GO
c.GO <- compareCluster(genelist, fun="enrichGO",OrgDb='org.Hs.eg.db')
dotplot(c.GO, showCategory=30)

# # GroupGO
# genedf.WT<-data.frame(Entrez=gene.WT,group="WT")
# genedf.KD<-data.frame(Entrez=gene.KD,group="KD")
# gene.df<-rbind(genedf.WT,genedf.KD)
# 
# c.GroupGO <- compareCluster(Entrez~group, data=gene.df,
#                              fun='groupGO', OrgDb='org.Hs.eg.db')
# 
# as.data.frame(c.GroupGO)
# dotplot(c.GroupGO, showCategory=30)
```

