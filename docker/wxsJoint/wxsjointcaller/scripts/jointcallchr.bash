#!/bin/bash
TMP_DIR=${OUTDIR}/tmp${LSB_JOBINDEX} && mkdir ${TMP_DIR}
DATABASE=$TMP_DIR/db
rm -R $DATABASE
${GATK4261} --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport \
    -L ${CHR} \
    --sample-name-map ${INDIR}/SampleMap.txt \
    --genomicsdb-workspace-path $DATABASE \
    --batch-size 50 \
    --reader-threads 6 \
    --genomicsdb-vcf-buffer-size 65536 \
    --tmp-dir ${TMP_DIR} \
&& ${GATK4261} --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs \
    -R ${REF_FASTA} \
    -V gendb://${DATABASE} \
    -O "${OUTDIR}/${LSB_JOBINDEX}.joint.vcf.gz" \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --tmp-dir ${TMP_DIR}
SUCCEEDED=$?
rm -R $TMP_DIR
exit $SUCCEEDED