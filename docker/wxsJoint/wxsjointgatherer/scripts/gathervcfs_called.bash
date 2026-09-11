#!/bin/bash
VCF_INPUTS=()
for ((int i; i < $NUM_INTERVALS; i++)); do
  VCF_INPUTS+="-I ${OUTDIR}/${INT}.joint.vcf.gz "
done

${GATK} \
  --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  GatherVcfs \
    -R ${REF_FASTA} \
    ${VCF_INPUTS[@]}\
    --TMP_DIR ${OUTDIR}/tmp \
    -RI \
    -O ${JOINT_VCF%.*.*}.unsorted.vcf.gz

${GATK} \
  --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SortVcf \
    -R ${REF_FASTA} \
    -I ${JOINT_VCF%.*.*}.unsorted.vcf.gz \
    --TMP_DIR ${OUTDIR}/tmp \
    -O ${JOINT_VCF}

rm ${JOINT_VCF%.*.*}.unsorted.vcf.gz
[[ -f ${JOINT_VCF}.tbi ]] || tabix $JOINT_VCF