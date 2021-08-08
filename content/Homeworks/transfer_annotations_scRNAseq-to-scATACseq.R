pbmc_atac <- 
pbmc_rna <- readRDS("data/PBMCs_10X/scRNAseq_pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
    reference = pbmc_rna,
    query = pbmc,
    reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
