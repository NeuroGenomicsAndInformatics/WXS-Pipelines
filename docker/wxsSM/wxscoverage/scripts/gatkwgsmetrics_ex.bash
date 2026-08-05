#!/bin/bash
CRM=$1
[[ -z $CRM ]] && CRM=${OUTDIR}/${CRAM}
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectWgsMetrics \
    -I $CRM \
    --INTERVALS ${REF_PADBED%.bed}.interval_list \
    -O ${CRM}.wgsmetrics_paddedexome.txt \
    -R ${REF_FASTA} \
    --TMP_DIR /tmp
