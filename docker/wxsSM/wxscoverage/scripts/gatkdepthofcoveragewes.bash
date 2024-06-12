#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  DepthOfCoverage \
    -I $1 \
    --INTERVALS ${REF_PADBED%.bed}.interval_list \
    -O $1.docmetrics_paddedexome \
    -R ${REF_FASTA} \
    --TMP_DIR /tmp \
    --omit-depth-output-at-each-base \
    --omit-interval-statistics
