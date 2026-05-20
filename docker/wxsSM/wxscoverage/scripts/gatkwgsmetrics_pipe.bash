#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectWgsMetrics \
    -I ${OUTDIR}/${CRAM} \
    -O ${OUTDIR}/${CRAM}.wgsmetrics.txt \
    -R ${REF_FASTA} \
    --TMP_DIR ${TMP_DIR} 
