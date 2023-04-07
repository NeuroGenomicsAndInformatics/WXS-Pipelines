#!/bin/bash
bash /scripts/splitintschr.bash
TMP_DIR=${OUTDIR}/tmp${LSB_JOBINDEX}
DATABASE=$TMP_DIR/db
INT_LISTS=($(ls /tmp | grep scattered))
${GATK4261} --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs \
    -R ${REF_FASTA} \
    -V gendb://${DATABASE} \
    -O "${OUTDIR}/${LSB_JOBINDEX}.joint.vcf.gz" \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir ${TMP_DIR}
SUCCEEDED=$?
rm -R $TMP_DIR
exit $SUCCEEDED