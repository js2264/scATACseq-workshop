tar -xzvf data/atac_pbmc_1k_v1_analysis.tar.gz
sed '1d' data/analysis/clustering/graphclust/clusters.csv | sed 's/,/\t/' > clusters.tsv
sinto filterbarcodes \
    -b data/PBMCs/atac_pbmc_1k_v1_possorted_bam.bam \
    --cells clusters.tsv \
    -p 12 

for FILE in *.bam
do

    SAMPLE=`echo "${FILE}" | sed 's,\..*,,'`
    FILTEREDBAM="${SAMPLE}"_filtered.bam
    TRACK="${SAMPLE}".bw

    samtools sort -@ 16 --output-fmt bam -n -T "${FILE}"_sorting "${FILE}" | \
    samtools fixmate -@ 16 --output-fmt bam -m - - | \
    samtools sort -@ 16 --output-fmt bam -T "${FILE}"_sorting2 - | \
    samtools view -@ 16 --output-fmt bam -f 2 -q 10 -1 -b - | \
    samtools sort -@ 16 --output-fmt bam -l 9 -T "${FILE}"_sorting3 -o "${FILTEREDBAM}"
    samtools index -@ 16 "${FILTEREDBAM}"

    Rscript -e "
        bam <- '${FILTEREDBAM}'
        track <- '${TRACK}'
        sample <- gsub('_filtered.*', '', bam)
        ncells <- sum(read.table('clusters.tsv')$V2 == sample)
        g <- VplotR::importPEBamFiles(
            bam, 
            genome = 'hg38', 
            shift_ATAC_fragments = TRUE, 
            verbose = FALSE
        )
        gr.cov <- IRanges::coverage(g) / (ncells*2)
        rtracklayer::export.bw(gr.cov, track)
    "

done