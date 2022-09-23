#!/bin/bash
${GATK} \
  --java-options "-Xmx30g -XX:ParallelGCThreads=1" \
  CollectVariantCallingMetrics \
    -I ${FINAL_OUTDIR}/${GVCF} \
    -O ${FINAL_OUTDIR}/${GVCF##*/}.vcfmetrics \
    -R ${REF_FASTA} \
    --DBSNP ${REF_DBSNP} \
    --THREAD_COUNT 6 \
    --GVCF_INPUT true \
    --TMP_DIR ${TMP_DIR}
