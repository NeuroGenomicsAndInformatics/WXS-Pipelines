#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectWgsMetrics \
    -I ${OUTDIR}/${CRAM} \
    --INTERVALS ${REF_PADBED%.bed}.interval_list \
    -O ${OUTDIR}/${CRAM}.wgsmetrics_paddedexome.txt \
    -R ${REF_FASTA} \
    --TMP_DIR ${TMP_DIR}
