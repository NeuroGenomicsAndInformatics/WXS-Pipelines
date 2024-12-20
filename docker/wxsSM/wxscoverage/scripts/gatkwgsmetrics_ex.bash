#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectWgsMetrics \
    -I $1 \
    --INTERVALS ${REF_PADBED%.bed}.interval_list \
    -O $1.wgsmetrics_paddedexome.txt \
    -R ${REF_FASTA} \
    --TMP_DIR /tmp
