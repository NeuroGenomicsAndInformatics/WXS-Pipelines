#!/bin/bash
CRM=$1
[[ -z $CRM ]] && CRM=${OUTDIR}/${CRAM}
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  DepthOfCoverage \
    -I $CRM \
    -L ${REF_PADBED%.bed}.interval_list \
    -O ${CRM}.docmetrics_paddedexome \
    -R ${REF_FASTA} \
    --omit-depth-output-at-each-base \
    --omit-interval-statistics
