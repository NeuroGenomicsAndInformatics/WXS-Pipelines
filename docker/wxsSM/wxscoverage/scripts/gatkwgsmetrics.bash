#!/bin/bash
${GATK} \
  --java-options "-Xmx70g -XX:ParallelGCThreads=1" \
  CollectWgsMetrics \
    -I ${FINAL_OUTDIR}/${CRAM} \
    -O /tmp/wgsmetrics.txt \
    -R ${REF_FASTA} \
    --TMP_DIR ${TMP_DIR} \
&& rsync /tmp/wgsmetrics.txt ${FINAL_OUTDIR}/${CRAM}.wgsmetrics.txt
