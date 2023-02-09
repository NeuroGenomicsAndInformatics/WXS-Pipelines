#!/bin/bash
bash /scripts/splitintschr.bash
TMP_DIR=${OUTDIR}/tmp${LSB_JOBINDEX} && mkdir ${TMP_DIR}
DATABASE=$TMP_DIR/db
rm -R $DATABASE
INT_LISTS=($(ls /tmp | grep scattered))
$GATK --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport \
    -L /tmp/${INT_LISTS[${LSB_JOBINDEX}-1]} \
    --sample-name-map ${INDIR}/SampleMap.txt \
    --genomicsdb-workspace-path $DATABASE \
    --batch-size 50 \
    --reader-threads 6 \
    --genomicsdb-vcf-buffer-size 65536 \
    --tmp-dir ${TMP_DIR} \
    --genomicsdb-shared-posixfs-optimizations true \
&& $GATK --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=2 -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs \
    -R ${REF_FASTA} \
    -V gendb://${DATABASE} \
    -O "${OUTDIR}/${LSB_JOBINDEX}.joint.vcf.gz" \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --tmp-dir ${TMP_DIR} \
    --genomicsdb-shared-posixfs-optimizations true 
SUCCEEDED=$?
rm -R $TMP_DIR
exit $SUCCEEDED