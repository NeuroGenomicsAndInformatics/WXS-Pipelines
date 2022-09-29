#!/bin/bash
${GATK} \
  --java-options "-Xmx70g -XX:ParallelGCThreads=1" \
  CollectWgsMetrics \
    -I ${FINAL_OUTDIR}/${CRAM} \
    -L ${REF_PADBED} \
    -O /tmp/wgsmetrics_paddedexome.txt \
    -R ${REF_FASTA} \
    --TMP_DIR ${TMP_DIR} \
&& rsync /tmp/wgsmetrics_paddedexome.txt ${FINAL_OUTDIR}/${CRAM}.wgsmetrics_paddedexome.txt
