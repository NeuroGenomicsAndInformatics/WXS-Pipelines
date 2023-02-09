#!/bin/bash
$GATK SplitIntervals \
    -R $REF_FASTA \
    -L /intlists/${CHR}.interval_list \
    -O /tmp \
    --scatter-count 50 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --subdivision-mode INTERVAL_COUNT
