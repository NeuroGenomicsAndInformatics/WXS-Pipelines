#!/bin/bash
for VAR in $(printenv | grep CUDA_VISIBLE_DEVICES); do
export ${VAR/CUDA/NVIDIA}
done
pbrun haplotypecaller \
  --ref ${REF_FASTA} \
  --in-bam ${OUTDIR}/${CRAM} \
  --in-recal-file ${OUTDIR}/${FULLSMID}.recal.txt \
  --out-variants ${OUTDIR}/${GVCF} \
  --gvcf \
  --num-gpus 1 \
  --tmp-dir ${TMP_DIR} \
  --annotation-group StandardAnnotation \
  --annotation-group StandardHCAnnotation \
  --annotation-group AS_StandardAnnotation \
&& rm -R $STAGE_INDIR