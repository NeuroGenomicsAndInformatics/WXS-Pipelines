#!/bin/bash
rsync -rL $STAGE_INDIR/ $INDIR
for BM in $(find $INDIR -name "*.bam"); do
  ln -s $BM /tmp/working.bam
  ${GATK} --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    RevertSamSpark \
    -I /tmp/working.bam \
    -O $TMP_DIR/reverted$LSB_JOBID.bam \
    --tmp-dir $TMP_DIR \
    --sort-order queryname \
    --read-validation-stringency SILENT \
    -- \
    --spark-master local[14]
  if [[ $? -ne 0 ]]; then
    rm -R $INDIR
    rm -R $TMP_DIR/reverted$LSB_JOBID.bam*
    exit 60
  fi
  rm $BM
  rm /tmp/working.bam
  mkdir $INDIR/${BM##*/}
  ${GATK} --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    SamToFastq \
    -I $TMP_DIR/reverted$LSB_JOBID.bam \
    --COMPRESS_OUTPUTS_PER_RG true \
    --OUTPUT_PER_RG true \
    --OUTPUT_DIR $INDIR/${BM##*/} \
    -RG_TAG PU \
    --TMP_DIR $TMP_DIR \
    --VALIDATION_STRINGENCY SILENT \
    --MAX_RECORDS_IN_RAM 10000000
  if [[ $? -ne 0 ]]; then
    rm -R $TMP_DIR/reverted$LSB_JOBID.bam*
    rm -R $INDIR
    exit 61
  fi
  rm $TMP_DIR/reverted$LSB_JOBID.bam
done
INFQ_FILE=${INDIR}/infqfile.txt
echo -n "" >$INFQ_FILE
for FQ in $(find $INDIR -name "*_1.f*q.gz"); do
  SM=$(echo $FULLSMID | cut -d^ -f1)
  BARCODE=$(echo $FULLSMID | cut -d^ -f2)
  PROJECT=$(echo $FULLSMID | cut -d^ -f3)
  FLOWCELL=$(zcat ${FQ} | head -n1 | cut -d: -f3)
  LANE=$(zcat ${FQ} | head -n1 | cut -d: -f4)
  FLOWLANE="${FLOWCELL}.${LANE}"
  [[ -z $FLOWCELL ]] && FLOWLANE=$(echo ${FQ##*/} | rev | cut -d_ -f2- | rev)
  echo "@RG\tID:${FLOWLANE}\tPL:illumina\tPU:${FLOWLANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >${OUTDIR}/${FULLSMID}.${FLOWLANE}.rgfile
  echo "${FQ} ${FQ/_1.f/_2.f} @RG\tID:${FLOWLANE}\tPL:illumina\tPU:${FLOWLANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >>${INFQ_FILE}
done
