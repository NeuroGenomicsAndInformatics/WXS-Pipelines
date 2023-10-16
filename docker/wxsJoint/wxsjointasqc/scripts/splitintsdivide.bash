#!/bin/bash
$GATK SplitIntervals \
    -R $REF_FASTA \
    -L /scratch1/fs1/cruchagac/WXSref/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    -O /tmp \
    --scatter-count $NUM_INTERVALS \
    --interval-merging-rule OVERLAPPING_ONLY \
    --subdivision-mode INTERVAL_SUBDIVISION
