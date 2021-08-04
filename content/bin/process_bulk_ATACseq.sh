SAMPLE=OmniATAC
R1=SRR891269_1.fastq.gz
R2=SRR891269_2.fastq.gz
GENOME=~/genomes/hg38/hg38.fa
GENOME=/data/hg38/hg38.fa
SAMFILE="${SAMPLE}".sam
FILTEREDBAM="${SAMPLE}"_filtered.bam
TRACK="${SAMPLE}".bw
PEAKFOLDER=peaks/"${SAMPLE}"
GENOMESIZE=2600000000
CPU=16

## -- Download reads 
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_2.fastq.gz

## -- Trim reads and run QC
trim_galore \
    --paired \
    --fastqc \
    --cores 4 \
    "${R1}" "${R2}"

## -- Map reads 
bowtie2 \
    --threads "${CPU}" \
    -x "${GENOME}" \
    --maxins 2000 \
    -1 `echo "${R1}" | sed 's,.fastq.gz,_val_1.fq.gz,'` \
    -2 `echo "${R2}" | sed 's,.fastq.gz,_val_2.fq.gz,'` > "${SAMFILE}"

## -- Filter reads 
samtools fixmate -@ "${CPU}" --output-fmt bam -m "${SAMFILE}" - | \
samtools sort -@ "${CPU}" --output-fmt bam -T "${SAMFILE}"_sorting - | \
samtools markdup -@ "${CPU}" --output-fmt bam -r -T "${SAMFILE}"_markdup - - | \
samtools view -@ "${CPU}" --output-fmt bam -f 2 -q 10 -1 -b - | \
samtools sort -@ "${CPU}" --output-fmt bam -l 9 -T "${SAMFILE}"_sorting2 -o "${FILTEREDBAM}"
samtools index -@ "${CPU}" "${FILTEREDBAM}"

## -- Generate track
bamCoverage \
    --bam "${FILTEREDBAM}" \
    --outFileName "${TRACK}" \
    --binSize 1 \
    --numberOfProcessors "${CPU}" \
    --normalizeUsing CPM \
    --skipNonCoveredRegions \
    --extendReads \
    --ignoreDuplicates

## -- Call peaks
macs2 callpeak \
    -t "${FILTEREDBAM}" \
    --format BAMPE \
    --gsize "${GENOMESIZE}" \
    --outdir "${PEAKFOLDER}" \
    --name "${SAMPLE}"
