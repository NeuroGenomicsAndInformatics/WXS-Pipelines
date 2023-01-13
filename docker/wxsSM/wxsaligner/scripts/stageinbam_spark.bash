#!/bin/bash
rsync -rL $STAGE_INDIR/ $INDIR
for BM in $(find $INDIR -name "*.bam"); do
ln -s $BM /tmp/working.bam
samtools index -@ $LSB_MAX_NUM_PROCESSORS /tmp/working.bam
${GATK} --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  RevertSamSpark \
    -I /tmp/working.bam \
    -O /tmp/reverted.bam \
    --tmp-dir $TMP_DIR \
    --sort-order queryname \
    --read-validation-stringency SILENT \
    -- \
    --spark-master local[14] \
&& rm $(find $INDIR -name "*.bam") \
&& ${GATK} --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SamToFastq \
    -I /tmp/reverted.bam \
    --COMPRESS_OUTPUTS_PER_RG true \
    --OUTPUT_PER_RG true \
    --OUTPUT_DIR $INDIR \
    -RG_TAG PU \
    --TMP_DIR $TMP_DIR \
    --VALIDATION_STRINGENCY SILENT \
    --MAX_RECORDS_IN_RAM 10000000
rm /tmp/*
done
sleep 10
INFQ_FILE=${INDIR}/infqfile.txt
echo -n "" > $INFQ_FILE
for FQ in $(find $INDIR -name "*_1.fastq.gz"); do
SM=$(echo $FULLSMID | cut -d^ -f1)
BARCODE=$(echo $FULLSMID | cut -d^ -f2)
PROJECT=$(echo $FULLSMID | cut -d^ -f3)
FLOWCELL=$(echo ${FQ##*/} | cut -d_ -f1 | cut -d. -f1)
LANE=$(echo ${FQ##*/} | cut -d_ -f1 | cut -d. -f2)
echo "@RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" > ${OUTDIR}/${FULLSMID}.${FLOWCELL}_${LANE}.rgfile
echo "${FQ} ${FQ/_1.fastq/_2.fastq} @RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >> ${INFQ_FILE}
done
