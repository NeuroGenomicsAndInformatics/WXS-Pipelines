#!/bin/bash
$GATK SplitIntervals \
    -R $REF_FASTA \
    -L /intlists/${CHR}.interval_list \
    -O /tmp \
    --scatter-count 100 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --subdivision-mode INTERVAL_SUBDIVISION
