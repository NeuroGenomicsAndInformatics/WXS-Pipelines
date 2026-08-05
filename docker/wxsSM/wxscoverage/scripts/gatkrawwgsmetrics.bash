#!/bin/bash
CRM=$1
[[ -z $CRM ]] && CRM=${OUTDIR}/${CRAM}
${GATK4261mod} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  CollectRawWgsMetrics \
    -I $CRM \
    -O ${CRM}.rawwgsmetrics.txt \
    -R ${REF_FASTA} \
    --TMP_DIR /tmp
