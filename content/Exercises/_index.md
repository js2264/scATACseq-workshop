---
title: "2. Exercises"
---

## Day 1

### Exercises session 1

- Check fragment size distribution in bulk ATAC-seq

```{r }
library(tidyverse)
library(plyranges)
fragments <- VplotR::importPEBamFiles(
    '?...?', 
    shift_ATAC_fragments = TRUE
)
df<-?...?(fragments) %>% 
    group_by(?...?) %>% 
    count() %>% 
    ungroup() 
p1<-ggplot(df, aes(x = ?...?, y = ?...?)) + 
    ?...?() + 
    labs(title = '?...?', x = '?...?', y = '?...?')
```

### Exercises session 2

- Check coverage of each cluster in `Loupe` at relevant loci to propose an annotation for each cluster. 

```shell
# Check e.g. `CD19`, `CD8`, `CD4`, `CD3`
```

## Day 2

### Exercises session 3

- Inspect genes closest to the most DA peaks (in DA analysis between clusters 10+16 and 4+5). Propose a function for cells in these groups.

```{r }
da_peaks <- readRDS('data/MouseBrain/da_peaks_10-16_vs_4-5.rds')
peaks <- rownames(da_peaks)[da_peaks$?...? <= 0.01 & da_peaks$?...? >= 1.5]
peaks <- GenomicRanges::GRanges(str_replace(peaks, '-', ':'))
Annotation(brain)[nearest(?...?, Annotation(?...?))]$?...?
```

### Exercises session 4

- Compare clustering to provided cell type activity for ce11 promoters. Comment. 

```{r }
dds <- readRDS('data/ATAC_worm/quantif.rds')
rlogs <- readRDS('?...?')
clusts <- cluster::pam(select(rlogs, -?...?), k = 9)
table(rowData(dds)$which.tissues[!rowData(dds)$allZero], clusts$?...?)
```

## Day 3

### Exercises session 5



### Exercises session 6

- Try different thresholds for tSNE embeddings. Comment. 

```{r }
set.seed(2021)
dev <- readRDS('data/scATAC_LSC-LMPP-mono/dev.rds')
p <- cowplot::plot_grid(
    chromVAR::plotDeviationsTsne(dev, chromVAR::deviationsTsne(dev, threshold = 1, perplexity = 20), annotation_name = "FOS", shiny = FALSE)[[1]],
    chromVAR::plotDeviationsTsne(dev, chromVAR::deviationsTsne(dev, threshold = 1.25, perplexity = 20), annotation_name = "FOS", shiny = FALSE)[[1]],
    chromVAR::plotDeviationsTsne(dev, chromVAR::deviationsTsne(dev, threshold = 1.5, perplexity = 20), annotation_name = "FOS", shiny = FALSE)[[1]],
    chromVAR::plotDeviationsTsne(dev, chromVAR::deviationsTsne(dev, threshold = 1.75, perplexity = 20), annotation_name = "FOS", shiny = FALSE)[[1]]
)
```
