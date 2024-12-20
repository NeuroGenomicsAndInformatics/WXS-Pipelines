#!/bin/bash
VCF_SHARDS=($(find ${OUTDIR} -maxdepth 1 -name "*joint.vcf.gz"))
echo ${#VCF_SHARDS[@]}
VCF_INPUTS=()
for VCF in ${VCF_SHARDS[@]}; do
  VCF_INPUTS+="-I ${VCF} "
done
${GATK} --java-options "-Xms200g -Xmx200g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SortVcf \
    -R ${REF_FASTA} \
    ${VCF_INPUTS[@]}\
    --TMP_DIR ${OUTDIR}/tmp \
    -O ${JOINT_VCF} \
&& rm -R $OUTDIR

for CHR in chr{1..22} chrX chrY; do
  bcftools view -r $CHR $JOINT_VCF > ${JOINT_VCF%.*.*.*.*}.$CHR.wgs.joint.vcf.gz
done