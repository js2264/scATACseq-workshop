---
title: "4. Prerequisites and required machine configuration"
---

## Prerequisites

The course is intended for those who have basic familiarity with `Unix` and the `R` scripting language.

* If a refresher is needed for Unix command line (hopefully not), please go over [this tutorial](https://ryanstutorials.net/linuxtutorial/) and its [companion cheatsheet](https://ryanstutorials.net/linuxtutorial/cheatsheet.php).
* Getting down to basics: an introduction to the fundamentals of R ([courtesy of Mark Ravinet](markravinet.github.io/Introduction.html)).
* Gentle introduction to `R/Biocondutor`: [here](https://bioconductor.github.io/BiocWorkshops/introduction-to-bioconductor-annotation-resources.html)
* For a full in-depth guide of `Bioconductor` ecosystem: read the comprehensive `R/Bioconductor` book from Kasper D. Hansen available under the CC BY-NC-SA 4.0 license [[PDF]](/{{<myPackageUrl>}}docs/bioconductor.pdf)

## Local configuration 

- Install R and other important softwares 

```sh
conda create --name scatac
conda activate scatac
conda install -c conda-forge -c bioconda \
    r-base r-essentials r-devtools r-tidyverse r-biocmanager r-polyclip r-rcpparmadillo \
    trim-galore bowtie2 samtools deeptools meme macs2 igv sinto yapc htseq
conda install -c anaconda hdf5
```

- The following R packages have also been installed (along with their many dependencies, of course!): 

```sh
Rscript -e "
## CRAN packages
install.packages('vroom', repos = 'https://cran.irsn.fr')
install.packages('umap', repos = 'https://cran.irsn.fr')
install.packages('corrplot', repos = 'https://cran.irsn.fr')
install.packages('gam', repos = 'https://cran.irsn.fr')
install.packages('ggbeeswarm', repos = 'https://cran.irsn.fr')
install.packages('ggthemes', repos = 'https://cran.irsn.fr')
install.packages('Matrix', repos = 'https://cran.irsn.fr')
install.packages('rgl', dependencies=TRUE, repos = 'https://cran.irsn.fr')
install.packages('hdf5r', dependencies=TRUE, repos = 'https://cran.irsn.fr')
install.packages('Signac', repos = 'https://cran.irsn.fr')
install.packages('Seurat', repos = 'https://cran.irsn.fr')
install.packages('BiocManager', repos = 'https://cran.irsn.fr')

## Bioconductor Packages
BiocManager::install('SingleCellExperiment', update = FALSE)
BiocManager::install('scran', update = FALSE)
BiocManager::install('scater', update = FALSE)
BiocManager::install('chromVAR', update = FALSE)
BiocManager::install('AnnotationHub', update = FALSE)
BiocManager::install('plyranges', update = FALSE)
BiocManager::install('VplotR', update = FALSE)
BiocManager::install('annotatr', update = FALSE)
BiocManager::install('rtracklayer', update = FALSE)
BiocManager::install('pheatmap', update = FALSE)
BiocManager::install('biovizBase', update = FALSE)
BiocManager::install('gprofiler2', update = FALSE)
BiocManager::install('DOSE', update = FALSE)
BiocManager::install('clusterProfiler', update = FALSE)
BiocManager::install('EnsDb.Mmusculus.v79', update = FALSE)
BiocManager::install('EnsDb.Hsapiens.v79', update = FALSE)
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene', update = FALSE)
BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene', update = FALSE)
BiocManager::install('org.Hs.eg.db', update = FALSE)
BiocManager::install('org.Mm.eg.db', update = FALSE)
"
```

- Additional dependencies: 

    * cellranger-atac

