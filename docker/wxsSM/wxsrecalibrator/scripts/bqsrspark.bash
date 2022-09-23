#!/bin/bash
ln -s ${OUTDIR}/$CRAM /tmp/working.cram
ln -s ${OUTDIR}/$CRAM.crai /tmp/working.cram.crai
${GATK} \
  --java-options "-Xmx100g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  BaseRecalibratorSpark \
    -I /tmp/working.cram \
    -R ${REF_FASTA} \
    --known-sites ${REF_MILLS_GOLD} \
    --known-sites ${REF_DBSNP} \
    --known-sites ${REF_ONEKGP1} \
    -O "/tmp/recal.txt" \
    -- \
    --spark-master local[$LSB_MAX_NUM_PROCESSORS]
cp /tmp/recal.txt ${OUTDIR}/${FULLSMID}.recal.txt
