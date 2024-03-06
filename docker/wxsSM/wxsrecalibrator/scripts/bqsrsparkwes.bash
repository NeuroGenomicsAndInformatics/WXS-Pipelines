#!/bin/bash
# Build Recal Table
THREADS=$(( LSB_MAX_NUM_PROCESSORS * 2 ))
ln -s ${OUTDIR}/$CRAM /tmp/working.cram
ln -s ${OUTDIR}/$CRAM.crai /tmp/working.cram.crai
${GATK} \
  --java-options "-Xmx40g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  BaseRecalibratorSpark \
    -I /tmp/working.cram \
    -R ${REF_FASTA} \
    -L ${REF_PADBED} \
    --known-sites ${REF_MILLS_GOLD} \
    --known-sites ${REF_DBSNP} \
    --known-sites ${REF_ONEKGP1} \
    -O "/tmp/recal.txt" \
    -- \
    --spark-master local[$THREADS]
cp /tmp/recal.txt ${OUTDIR}/${FULLSMID}.recal.txt

# Apply Recal Table
${GATK} \
  --java-options "-Xmx40g -XX:ParallelGCThreads=1" \
  ApplyBQSR \
    -I ${OUTDIR}/${BAM} \
    -bqsr ${OUTDIR}/${FULLSMID}.recal.txt \
    -R ${REF_FASTA} \
    -L ${REF_PADBED} \
    -O ${OUTDIR}/${FULLSMID}.recal.bam
samtools index ${OUTDIR}/${FULLSMID}.recal.bam