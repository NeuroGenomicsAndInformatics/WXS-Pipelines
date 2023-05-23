#!/bin/bash
OUTDIR="$1"
JOINT_VCF="$2"
VCF_SHARDS=($(find ${OUTDIR} -maxdepth 1 -name "*joint.vcf.gz"))
echo ${#VCF_SHARDS[@]}
VCF_INPUTS=()
for VCF in ${VCF_SHARDS[@]}; do
  VCF_INPUTS+="-I ${VCF} "
done
${GATK} --java-options "-Xms50g -Xmx50g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SortVcf \
    -R /03-DryLab/01-References_Software/01-References/software_refs/GATK_pipeline_files/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    ${VCF_INPUTS[@]}\
    --TMP_DIR ${OUTDIR}/tmp \
    -O ${JOINT_VCF}
CODE=$?
rm -R $OUTDIR
exit $CODE