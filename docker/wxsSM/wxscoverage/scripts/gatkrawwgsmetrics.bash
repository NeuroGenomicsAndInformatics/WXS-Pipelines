#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectRawWgsMetrics \
    -I $1 \
    -O $1.rawwgsmetrics.txt \
    -R ${REF_FASTA} \
    --TMP_DIR /tmp
