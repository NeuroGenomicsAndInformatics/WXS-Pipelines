#!/bin/bash
$GATK SplitIntervals \
    -R $REF_FASTA \
    -L /intlists/WGS.interval_list \
    -O /tmp \
    --scatter-count $NUM_INTERVALS \
    --interval-merging-rule OVERLAPPING_ONLY \
    --subdivision-mode INTERVAL_SUBDIVISION
