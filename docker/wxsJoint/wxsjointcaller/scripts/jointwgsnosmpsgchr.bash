#!/bin/bash
COHORT=$1
OUTDIR=$2
INDIR=/storage1/fs1/cruchagac/Active/matthew.j/c1in/${COHORT}
mkdir ${OUTDIR}/tmp
TMP_DIR=${OUTDIR}/tmp/tmp${LSB_JOBINDEX} && mkdir ${TMP_DIR}
DATABASE=$TMP_DIR/db
rm -R $DATABASE
INT_LISTS=($(ls ${OUTDIR}/intlists))
/ref/gatk-4.2.6.1_mod/gatk --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport \
    -L ${OUTDIR}/intlists/${INT_LISTS[${LSB_JOBINDEX}-1]} \
    --sample-name-map ${INDIR}/SampleMap.txt \
    --genomicsdb-workspace-path $DATABASE \
    --batch-size 100 \
    --tmp-dir ${TMP_DIR} \
    --genomicsdb-shared-posixfs-optimizations true \
&& /ref/gatk-4.2.6.1_mod/gatk --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs \
    -R /ref/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    -V gendb://${DATABASE} \
    -O "${OUTDIR}/jointvcfs/${COHORT}.${OUTDIR##*/}.${LSB_JOBINDEX}.joint.g.vcf.gz" \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --tmp-dir ${TMP_DIR} \
    --genomicsdb-shared-posixfs-optimizations true \
&& rm -R ${TMP_DIR}
