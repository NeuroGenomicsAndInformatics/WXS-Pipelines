#!/bin/bash
rsync -rL $STAGE_INDIR/ $INDIR
pbrun bam2fq \
  --ref /ref/$1 \
  --in-bam $(find $INDIR -name "*.cram") \
  --out-prefix ${INDIR} \
  --rg-tag PU \
  --num-gpus 1 \
  --tmp-dir ${TMP_DIR}
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
rm $(find $INDIR -name "*.cram")
