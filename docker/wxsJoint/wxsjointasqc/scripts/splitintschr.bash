#!/bin/bash
CHR=$1
$GATK SplitIntervals \
    -R $REF_FASTA \
    -L /ref/intlists/${CHR}/${CHR}.interval_list \
    -O /tmp \
    --scatter-count 50 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --subdivision-mode INTERVAL_SUBDIVISION
