---
title: "2. Program"
---

## Day 1: Processing bulk and single-cell ATAC-seq reads

### ATAC-seq experimental approaches 
>30'

  - bulk ATAC-seq methods: classical, omniATAC-seq, low input ATAC-seq, ...
  - "single cell" ATAC-seq methods: droplet-based scATACseq, indexed scATACseq, joint scRNAseq and scATACseq, ...

### Processing of bulk ATAC-seq
>45'

  - Sequencing QC
  - Mapping reads and correcting
  - Visualizing and QC-ing the results

  ### Processing of single-cell ATAC-seq
>45'

  - 10X scATACseq: from reads to matrix counts with `cellranger`
  - QC-ing with `Loupe`

### Data output: loci of accessible chromatin 
>Extra

  - Distribution of fragment sizes over ATAC peaks
  - Distance to closest TSS
  - Annotating closest genes 

## Day 2: Peak-centered downstream analysis 

### Bulk ATACseq and peak sets
>45'

  - Integrating peak sets from different experiments
  - Quantifying accessibility at each ATAC peak
  - Differential accessibility tests
  - Clustering peak sets

### Clustering scATACseq cells
>45'

  - Reading `cellranger` output in `R` 
  - Normalizing scATACseq counts
  - Clustering cells by their accessible loci
  - Annotating cell types

### Insigths into regulatory functions 
>30' 

  - Interpreting chromatin context(s) at accessible loci
  - Interpreting DNA sequence context(s) at accessible loci
  - Leveraging Bioconductor `AnnotationHub` and `ExperimentHub` resources

### Insigths into biological functions
>Extra

  - Over-representation analyses
  - Gene Set Enrichment Analyses (GSEA)

## Day 3: "Multi-omics" integration 

### Overlapping scATACseq on top of scRNAseq cells
>30'

### Processing joint scRNAseq + scATACseq
>45'

### Building regulatory networks from scATACseq
>45'

  - Gene Regulatory Networks 
  - CICERO


