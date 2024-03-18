#!/bin/bash
rsync -rL $STAGE_INDIR/ $INDIR
for BM in $(find $INDIR -name "*.bam"); do
  mkdir ${BM%.bam}
  ${GATK} --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    RevertSam \
    -I $BM \
    -O /dev/stdout \
    --TMP_DIR $TMP_DIR \
    -SO queryname \
    --VALIDATION_STRINGENCY SILENT |
    ${GATK} --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
      SamToFastq \
      -I /dev/stdin \
      --COMPRESS_OUTPUTS_PER_RG true \
      --OUTPUT_PER_RG true \
      --OUTPUT_DIR ${BM%.bam} \
      -RG_TAG ID \
      --TMP_DIR $TMP_DIR \
      --VALIDATION_STRINGENCY SILENT \
      --MAX_RECORDS_IN_RAM 10000000
  if [[ $? -ne 0 ]]; then
    rm -R $INDIR
    exit 61
  fi
  rm $BM
done
INFQ_FILE=${INDIR}/infqfile.txt
echo -n "" > $INFQ_FILE
for FQ in $(find $INDIR -name "*1.f*q.gz"); do
SM=$(echo $FULLSMID | cut -d^ -f1)
BARCODE=$(echo $FULLSMID | cut -d^ -f2)
PROJECT=$(echo $FULLSMID | cut -d^ -f3)
HEAD=$(zcat ${FQ} | head -n1)
HEAD_CHECK="${HEAD//[^:]}"
if [[ ${#HEAD_CHECK} -eq 4 ]]; then
FLOWCELL=$(zcat ${FQ} | head -n1 | cut -d: -f1 | cut -d '@' -f2)
LANE=$(zcat ${FQ} | head -n1 | cut -d: -f2)
FLOWLANE="${FLOWCELL}.${LANE}"
elif [[ ${#HEAD_CHECK} -gt 4 ]]; then
FLOWCELL=$(zcat ${FQ} | head -n1 | cut -d: -f3)
LANE=$(zcat ${FQ} | head -n1 | cut -d: -f4)
FLOWLANE="${FLOWCELL}.${LANE}"
else
FLOWLANE=$(echo ${FQ##*/} | rev | cut -d_ -f2- | rev)
fi
echo "@RG\tID:${FLOWLANE}\tPL:illumina\tPU:${FLOWLANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" > ${OUTDIR}/${FULLSMID}.${FLOWLANE}.rgfile
echo "${FQ} ${FQ/_1.f/_2.f} @RG\tID:${FLOWLANE}\tPL:illumina\tPU:${FLOWLANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >> ${INFQ_FILE}
done