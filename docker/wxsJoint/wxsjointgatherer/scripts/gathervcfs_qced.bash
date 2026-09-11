#!/bin/bash
VCF_INPUTS=()
for (( i=0;i<$NUM_INTERVALS;i++)); do 
  VCF_INPUTS+="-I $(find $OUTDIR -name "$i.*ABfiltered.vcf.gz") "
done

${GATK} \
  --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  GatherVcfs \
    -R ${REF_FASTA} \
    ${VCF_INPUTS[@]}\
    --TMP_DIR ${OUTDIR}/tmp \
    -RI \
    -O ${JOINT_VCF%.*.*}.qced.unsorted.vcf.gz

${GATK} \
  --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SortVcf \
    -R ${REF_FASTA} \
    -I ${JOINT_VCF%.*.*}.qced.unsorted.vcf.gz \
    --TMP_DIR ${OUTDIR}/tmp \
    -O ${JOINT_VCF%.*.*}.qced.vcf.gz \
&& rm ${JOINT_VCF%.*.*}.qced.unsorted.vcf.gz \
&& find ${OUTDIR} -maxdepth 1 -name "*joint.vcf.gz*" -delete \
&& find ${OUTDIR} -maxdepth 1 -name "*.calc" -delete \
&& find ${OUTDIR} -maxdepth 1 -name "*.table" -delete \
&& find ${OUTDIR} -maxdepth 1 -name "*recalibrated.vcf.gz*" -delete

[[ -f ${JOINT_VCF%.*.*}.qced.vcf.gz.tbi ]] || tabix ${JOINT_VCF%.*.*}.qced.vcf.gz