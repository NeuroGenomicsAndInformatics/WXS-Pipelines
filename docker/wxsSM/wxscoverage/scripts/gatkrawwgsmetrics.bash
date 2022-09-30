#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectRawWgsMetrics \
    -I ${FINAL_OUTDIR}/${CRAM} \
    -O ${FINAL_OUTDIR}/${CRAM}.rawwgsmetrics.txt \
    -R ${REF_FASTA} \
    --TMP_DIR ${TMP_DIR}
