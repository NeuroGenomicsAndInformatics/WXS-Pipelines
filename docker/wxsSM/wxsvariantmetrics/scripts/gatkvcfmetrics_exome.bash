#!/bin/bash
VCF=$1
[[ -z $VCF ]] && VCF=${OUTDIR}/${GVCF}
${GATK} \
  --java-options "-Xmx10g -XX:ParallelGCThreads=1" \
  CollectVariantCallingMetrics \
    -I $VCF \
    -O ${VCF}.vcfmetrics \
    -R ${REF_FASTA} \
    --DBSNP ${REF_DBSNP} \
    --THREAD_COUNT 3 \
    --TARGET_INTERVALS ${REF_PADBED%.bed}.interval_list \
    --GVCF_INPUT true \
    --TMP_DIR /tmp
