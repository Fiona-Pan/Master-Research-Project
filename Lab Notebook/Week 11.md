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

