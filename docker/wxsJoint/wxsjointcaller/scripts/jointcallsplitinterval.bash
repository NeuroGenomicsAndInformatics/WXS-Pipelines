#!/bin/bash
INTERVAL=$1
bash /scripts/splitintsdivide.bash
TMP_DIR=${OUTDIR}/tmp${INTERVAL} && mkdir ${TMP_DIR}
DATABASE=$TMP_DIR/db
rm -R $DATABASE
INT_LISTS=($(ls /tmp | grep scattered))
${GATK4261} --java-options "-Xms100g -Xmx100g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport \
    -L /tmp/${INT_LISTS[$INTERVAL]} \
    --sample-name-map ${INDIR}/SampleMap.txt \
    --genomicsdb-workspace-path $DATABASE \
    --batch-size 50 \
    --consolidate true \
    --genomicsdb-shared-posixfs-optimizations true \
    --genomicsdb-vcf-buffer-size 65536 \
    --tmp-dir ${TMP_DIR} \
&& ${GATK4261} --java-options "-Xms100g -Xmx100g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs \
    -R ${REF_FASTA} \
    -V gendb://${DATABASE} \
    -O "${OUTDIR}/${INTERVAL}.joint.vcf.gz" \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir ${TMP_DIR}
SUCCEEDED=$?
rm -R $TMP_DIR
exit $SUCCEEDED