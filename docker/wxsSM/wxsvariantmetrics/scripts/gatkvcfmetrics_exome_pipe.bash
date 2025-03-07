#!/bin/bash
${GATK} \
  --java-options "-Xmx10g -XX:ParallelGCThreads=1" \
  CollectVariantCallingMetrics \
    -I ${FINAL_OUTDIR}/${GVCF} \
    -O ${FINAL_OUTDIR}/${GVCF##*/}.vcfmetrics \
    -R ${REF_FASTA} \
    --DBSNP ${REF_DBSNP} \
    --TARGET_INTERVALS ${REF_PADBED%.bed}.interval_list \
    --THREAD_COUNT 3 \
    --GVCF_INPUT true \
    --TMP_DIR ${TMP_DIR}
