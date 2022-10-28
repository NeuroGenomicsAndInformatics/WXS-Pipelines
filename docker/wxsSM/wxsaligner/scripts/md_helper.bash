#!/bin/bash
MD_INPUTS=()
for BM in $(find $INDIR -name "*.fastq*.bam"); do
MD_INPUTS+="-I ${BM} "
done
${GATK} \
  --java-options "-Xmx170g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  MarkDuplicates \
    ${MD_INPUTS[@]}\
    -O ${OUTDIR}/${CRAM} \
    -M ${METDIR}/${FULLSMID}.dup.metrics.txt \
    -R ${REF_FASTA} \
    --SORTING_COLLECTION_SIZE_RATIO 0.25 \
    --TMP_DIR ${TMP_DIR} \
    --MAX_RECORDS_IN_RAM 2500000
rm -R ${INDIR}
samtools index -@ $LSB_MAX_NUM_PROCESSORS $OUTDIR/$CRAM
