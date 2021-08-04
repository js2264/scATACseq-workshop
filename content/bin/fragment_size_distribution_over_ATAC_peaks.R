fragments <- importPEBamFiles(
    'path/to/ATAC.bam', 
    shift_ATAC_fragments = TRUE
)
peaks <- import('path/to/peaks.narrowPeak')

## -- Fragment size distribution
df<-subsetByOverlaps(fragments, peaks) %>% 
    as_tibble() %>% 
    group_by(width) %>% 
    count() %>% 
    ungroup() 
p1<-ggplot(df, aes(x = width, y = n)) + 
    geom_line() + 
    theme_minimal() + 
    theme(legend.position = 'none') + 
    labs(title = 'Distribution of ATAC-seq fragment size', x = 'Fragment size', y = '# of fragments within peaks')
