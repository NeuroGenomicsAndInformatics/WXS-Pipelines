#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  DepthOfCoverage \
    -I $1 \
    -L ${REF_PADBED%.bed}.interval_list \
    -O $1.docmetrics_paddedexome \
    -R ${REF_FASTA} \
    --omit-depth-output-at-each-base \
    --omit-interval-statistics
