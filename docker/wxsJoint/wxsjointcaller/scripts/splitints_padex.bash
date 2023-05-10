#!/bin/bash
$GATK SplitIntervals \
    -R $REF_FASTA \
    -L /ref/20190812_GATK_38_googlebundle/Capture_Padded.GRCh38.interval_list \
    -O /tmp \
    --scatter-count $1 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
