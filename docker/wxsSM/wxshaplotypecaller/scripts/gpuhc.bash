#!/bin/bash
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
  --annotation-group AS_StandardAnnotation
