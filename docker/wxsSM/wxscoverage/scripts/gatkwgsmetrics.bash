#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectWgsMetrics \
    -I $1 \
    -O $1.wgsmetrics.txt \
    -R ${REF_FASTA} \
    --TMP_DIR /tmp
