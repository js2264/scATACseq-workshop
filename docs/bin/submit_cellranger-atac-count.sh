#!/usr/bin/env bash

module load cellranger-atac

cellranger-atac count \
    --id "${ID}" \
    --reference "${REF}" \
    --fastqs "${FASTQS}" \
    --localcores "${CPU}" \
    --localmem "${MEM}"

