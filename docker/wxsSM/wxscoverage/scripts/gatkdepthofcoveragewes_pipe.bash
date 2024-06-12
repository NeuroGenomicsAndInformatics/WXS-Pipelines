#!/bin/bash
${GATK} \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  DepthOfCoverage \
    -I ${FINAL_OUTDIR}/${CRAM} \
    -L ${REF_PADBED%.bed}.interval_list \
    -O ${FINAL_OUTDIR}/${CRAM}.docmetrics_paddedexome \
    -R ${REF_FASTA} \
    --omit-depth-output-at-each-base \
    --omit-interval-statistics
