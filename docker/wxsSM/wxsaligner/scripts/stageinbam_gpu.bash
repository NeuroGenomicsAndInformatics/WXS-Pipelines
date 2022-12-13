#!/bin/bash
for VAR in $(printenv | grep CUDA_VISIBLE_DEVICES); do
export ${VAR/CUDA/NVIDIA}
done
rsync -rL $STAGE_INDIR/ $INDIR
for BM in $(find $INDIR -name "*.bam"); do
samtools index -@ $LSB_MAX_NUM_PROCESSORS $BM
pbrun bam2fq \
  --in-bam $(find $INDIR -name "*.bam") \
  --out-prefix ${INDIR}/${FULLSMID} \
  --rg-tag PU \
  --tmp-dir ${TMP_DIR}
done
INFQ_FILE=${INDIR}/infqfile.txt
echo -n "" > $INFQ_FILE
for FQ in $(find $INDIR -name "*_1.fastq*"); do
SM=$(echo $FULLSMID | cut -d^ -f1)
BARCODE=$(echo $FULLSMID | cut -d^ -f2)
PROJECT=$(echo $FULLSMID | cut -d^ -f3)
FLOWCELL=$(echo ${FQ##*/} | cut -d_ -f2 | cut -d. -f1)
LANE=$(echo ${FQ##*/} | cut -d_ -f2 | cut -d. -f2)
echo "@RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" > ${OUTDIR}/${FULLSMID}.${FLOWCELL}_${LANE}.rgfile
echo "${FQ} ${FQ/_1.fastq/_2.fastq} @RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >> ${INFQ_FILE}
done
rm $(find $INDIR -name "*.cram")
