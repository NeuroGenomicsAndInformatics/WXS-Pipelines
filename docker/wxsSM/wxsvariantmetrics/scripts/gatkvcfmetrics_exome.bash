#!/bin/bash
if [ -z $FULLSMID ]; then
${GATK} \
  --java-options "-Xmx10g -XX:ParallelGCThreads=1" \
  CollectVariantCallingMetrics \
    -I $1 \
    -O $1.vcfmetrics \
    -R ${REF_FASTA} \
    --DBSNP ${REF_DBSNP} \
    --THREAD_COUNT 3 \
    --TARGET_INTERVALS ${REF_PADBED%.bed}.interval_list \
    --GVCF_INPUT true \
    --TMP_DIR /tmp
else
  bash /scripts/gatkvcfmetrics_exome_pipe.bash
fi
