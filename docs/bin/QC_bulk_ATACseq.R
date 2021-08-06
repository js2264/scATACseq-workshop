library(AnnotationDbi)
library(tidyverse)
library(plyranges)
library(annotatr)
library(VplotR)

## -- Importing genomic features
hg_genes <- AnnotationHub::AnnotationHub()[['AH92109']] %>% 
    filter(type == 'gene', gene_biotype == 'protein_coding') %>% 
    filter(!(grepl('GL|KI', seqnames)))
seqlevelsStyle(hg_genes) <- 'UCSC'
hg_TSSs <- hg_genes %>% 
    anchor_start() %>% 
    resize(width = 1) %>% 
    unanchor()
hg_promoters <- hg_TSSs %>%
    flank_upstream(4000) %>% 
    shift_downstream(2000) %>% 
    mutate(ID = paste0('prom_', 1:n()))
hg_features <- build_annotations(
    genome = 'hg38',
    annotations = 'hg38_basicgenes'
)
hg_chromstates <- AnnotationHub::AnnotationHub()[['AH46969']]

## -- Importing fragments and peaks as GRanges
fragments <- importPEBamFiles(
    'OmniATAC_filtered.bam', 
    shift_ATAC_fragments = TRUE
)
peaks <- import('peaks/OmniATAC/OmniATAC_peaks.narrowPeak')
peak_summits <- import('peaks/OmniATAC/OmniATAC_summits.bed')

## -- Fragment size distribution
df<-as_tibble(fragments) %>% 
    group_by(width) %>% 
    count() %>% 
    ungroup() 
p1<-ggplot(df, aes(x = width, y = n)) + 
    geom_line() + 
    theme_minimal() + 
    theme(legend.position = 'none') + 
    labs(title = 'Distribution of ATAC-seq fragment size', x = 'Fragment size', y = '# of fragments')

## -- Vplot
p2<-plotVmat(
    fragments, 
    hg_TSSs, 
    normFun = 'zscore', 
    ylim = c(20, 800), 
    xlim = c(-500, 500)
)

## -- Distance to TSS
df<-add_nearest_distance(fragments, hg_TSSs) %>%
    as_tibble() %>% 
    mutate(distance = cut(distance, breaks = seq(0, 5e5, by = 10), include.lowest = TRUE) %>% as.numeric() %>% `*`(10)) %>%
    group_by(distance) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(cumsum = cumsum(n), pct = cumsum/max(cumsum))
p3<-ggplot(df, aes(x = distance, y = pct)) + 
    geom_line() + 
    theme_minimal() + 
    scale_x_log10() +
    theme(legend.position = 'none') + 
    labs(title = 'Cumulative distribution of ATAC-seq fragments', x = 'Distance from TSS', y = 'Cum. %')

## -- Promoter (TSS) Enrichment Score
df<-resize(fragments, width = 1, fix = 'center') %>% 
    join_overlap_left(hg_promoters) %>% 
    plyranges::select(ID) %>%
    as_tibble() %>% 
    filter(!is.na(ID)) %>%
    left_join(as_tibble(hg_promoters) %>% dplyr::select(start, end, ID), by = 'ID') %>% 
    mutate(
        prom_mid = start.y + (end.y - start.y) / 2, 
        distance = start.x - prom_mid, 
        binned_distance = round(distance/100, 0)*100
    ) %>%
    count(binned_distance) %>%
    mutate(
        bg = mean(.data$n[abs(.data$binned_distance) == 2000]), 
        enrich_score = n/bg
    )
TSSES <- round(df[df$binned_distance == 0, 'enrich_score'], 2)
p4<-ggplot(df, aes(x = binned_distance, y = enrich_score)) + 
    geom_line() + 
    theme_minimal() + 
    theme(legend.position = 'none') + 
    labs(title = glue::glue('TSS enrichment score (Signal-to-noise ratio): {TSSES}'), x = 'Distance from TSS', y = 'Enrich. score')

## -- Location of peaks vs genomic features
peaks_annotated <- annotate_regions(
    peak_summits, 
    annotations = hg_features,
    ignore.strand = TRUE
)
p5<-plot_annotation(
    annotated_regions = peaks_annotated,
    annotation_order = paste0('hg38_genes_', c('promoters', '1to5kb', '5UTRs', 'exons', 'introns', '3UTRs')),
    plot_title = 'ATAC peaks',
    x_label = 'Genomic features',
    y_label = 'Count'
) + theme_minimal()

## -- FRiP 
df<-add_nearest_distance(fragments, peaks) %>%
    as_tibble() %>% 
    mutate(isInPeak = distance == 0)
pct <- round(sum(df$isInPeak == 1) / nrow(df) * 100, 2)
p6<-ggplot(df, aes(y = isInPeak)) + 
    geom_bar() + 
    theme_minimal() + 
    theme(legend.position = 'none') + 
    labs(title = glue::glue('Fraction of reads in peaks (FRiP): {pct}'), x = '# of fragments', y = 'Fragment within peak')

## -- All plots together
p <- cowplot::plot_grid(
    p1, p2, p3, p4, p5, p6, 
    nrow = 2
)
ggsave('ATACseq_QC.pdf', w = 15, h = 10)