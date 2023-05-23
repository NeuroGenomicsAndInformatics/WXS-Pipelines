#!/bin/bash
INDIR=$1
OUTDIR=$2
INTERVAL=$3

JOBID=$(echo ${INTERVAL##*/} | cut -d. -f1)

REF_DIR=/03-DryLab/01-References_Software/01-References/software_refs/GATK_pipeline_files
TMP_DIR=${OUTDIR}/tmp$JOBID && mkdir ${TMP_DIR}
DATABASE=$TMP_DIR/db
rm -R $DATABASE
gatk --java-options "-Xms40g -Xmx40g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenomicsDBImport \
    -L $INTERVAL \
    --sample-name-map ${INDIR}/SampleMap.txt \
    --genomicsdb-workspace-path $DATABASE \
    --batch-size 50 \
    --consolidate true \
    --tmp-dir ${TMP_DIR} \
&& gatk --java-options "-Xms40g -Xmx40g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    GenotypeGVCFs \
    -R ${REF_DIR}/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    -V gendb://${DATABASE} \
    -O "${OUTDIR}/${JOBID}.joint.vcf.gz" \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --tmp-dir ${TMP_DIR}

rm -R ${TMP_DIR}