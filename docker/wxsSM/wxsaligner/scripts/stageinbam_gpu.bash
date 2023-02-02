#!/bin/bash
for VAR in $(printenv | grep CUDA_VISIBLE_DEVICES); do
export ${VAR/CUDA/NVIDIA}
done
rsync -rL $STAGE_INDIR/ $INDIR
samtools merge -@ 7 -cu -l 0 $INDIR/merged.bam $INDIR/*.bam
for BM in $(find $INDIR -name "*.bam" | grep -v merged); do rm $BM; done
pbrun bam2fq \
  --in-bam $INDIR/merged.bam \
  --out-prefix ${INDIR}/${FULLSMID} \
  --rg-tag PU \
  --tmp-dir ${TMP_DIR}
rm $INDIR/*.bam
sleep 10
INFQ_FILE=${INDIR}/infqfile.txt
echo -n "" > $INFQ_FILE
for FQ in $(find $INDIR -name "*_1.fastq*"); do
SM=$(echo $FULLSMID | cut -d^ -f1)
BARCODE=$(echo $FULLSMID | cut -d^ -f2)
PROJECT=$(echo $FULLSMID | cut -d^ -f3)
FLOWCELL=$(echo ${FQ##*/} | rev | cut -d_ -f2- | rev | cut -d. -f1)
LANE=$(echo ${FQ##*/} | rev | cut -d_ -f2- | rev | cut -d. -f2)
echo "@RG\tID:${FLOWCELL}.${LANE}\tPL:illumina\tPU:${FLOWCELL}.${LANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" > ${OUTDIR}/${FULLSMID}.${FLOWCELL}_${LANE}.rgfile
echo "${FQ} ${FQ/_1.fastq/_2.fastq} @RG\tID:${FLOWCELL}.${LANE}\tPL:illumina\tPU:${FLOWCELL}.${LANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >> ${INFQ_FILE}
done
