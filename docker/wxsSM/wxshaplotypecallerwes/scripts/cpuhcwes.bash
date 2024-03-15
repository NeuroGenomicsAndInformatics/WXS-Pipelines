#!/bin/bash
${GATK} \
  --java-options "-Xmx40g -XX:ParallelGCThreads=1" \
  HaplotypeCaller \
    -I ${OUTDIR}/${FULLSMID}.recal.cram \
    -R ${REF_FASTA} \
    --dbsnp ${REF_DBSNP} \
    -L ${REF_PADBED} \
    -ERC GVCF \
    -O ${OUTDIR}/${GVCF} \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
&& rm -R ${STAGE_INDIR}
rm ${OUTDIR}/${FULLSMID}.recal.cram*