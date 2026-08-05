#!/bin/bash
VCF_SHARDS=($(find ${OUTDIR} -maxdepth 1 -name "*joint.vcf.gz"))
echo ${#VCF_SHARDS[@]}
VCF_INPUTS=()
for VCF in ${VCF_SHARDS[@]}; do
  VCF_INPUTS+="-I ${VCF} "
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