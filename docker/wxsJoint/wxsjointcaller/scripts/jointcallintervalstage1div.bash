#!/bin/bash
bash /scripts/splitintschrdivide.bash
TMP_DIR=${OUTDIR}/tmp${LSB_JOBINDEX} && mkdir ${TMP_DIR}
DATABASE=$TMP_DIR/db
rm -R $DATABASE
INT_LISTS=($(ls /tmp | grep scattered))
${GATK4261} --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport \
    -L /tmp/${INT_LISTS[${LSB_JOBINDEX}-1]} \
    --interval-padding 200 \
    --sample-name-map ${INDIR}/SampleMap.txt \
    --genomicsdb-workspace-path $DATABASE \
    --batch-size 150 \
    --consolidate true \
    --genomicsdb-shared-posixfs-optimizations true \
    --genomicsdb-vcf-buffer-size 655360 \
    --tmp-dir ${TMP_DIR}
