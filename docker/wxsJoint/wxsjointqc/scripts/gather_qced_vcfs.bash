#!/bin/bash
# This script gathers sites-only vcfs.
INPUT_DIR=$1

# Gather QCed vcfs
QC_INPUTS=()
for QCED_VCF in ${INPUT_DIR}/${INPUT_DIR##*/}.{1..50}.*.splitMA.geno.vcf.gz; do \
    QC_INPUTS+=" -I $QCED_VCF"
done
OUTPUT_QC_VCF="${INPUT_DIR}/${INPUT_DIR##*/}.qced.gathered.vcf.gz"

$GATK GatherVcfs \
    -R ${REF_FASTA} \
    --COMPRESSION_LEVEL 7 \
    ${QC_INPUTS[@]} \
    -O $OUTPUT_QC_VCF
tabix -f $OUTPUT_QC_VCF
SUCCESS_QC=$?

# Gather sites-only annotated vcfs
ANN_INPUTS=()
for ANNED_VCF in ${INPUT_DIR}/${INPUT_DIR##*/}.{1..50}.*.splitMA.*.ann.vcf.gz; do \
    ANN_INPUTS+=" -I $ANNED_VCF"
done
OUTPUT_ANN_VCF="${INPUT_DIR}/${INPUT_DIR##*/}.qced.gathered.ann.vcf.gz"

$GATK GatherVcfs \
    -R ${REF_FASTA} \
    --COMPRESSION_LEVEL 7 \
    ${ANN_INPUTS[@]} \
    -O $OUTPUT_ANN_VCF
tabix -f $OUTPUT_ANN_VCF
SUCCESS_ANN=$?

[[ $SUCCESS_QC -eq 0 ]] \
    && rm ${INPUT_DIR}/${INPUT_DIR##*/}.{1..50}.*.geno.vcf.gz* \

[[ $SUCCESS_ANN -eq 0 ]] \
    && rm ${INPUT_DIR}/${INPUT_DIR##*/}.{1..50}.*.ann.vcf.gz*

SUCCESS=66
[[ $SUCCESS_QC -eq 0 ]] && [[ $SUCCESS_ANN -eq 0 ]] \
    && rm ${INPUT_DIR}/${INPUT_DIR##*/}.{bed,interval_list,vcf.bgz,vcf.bgz.csi} \
    && rm -R ${INPUT_DIR}/intlists \
    && SUCCESS=0



# Gather count files
[[ -f ${INPUT_DIR}/${INPUT_DIR##*/}_combined_counts.csv ]] || /scripts/GatherCounts.R $INPUT_DIR

[[ -f ${INPUT_DIR}/${INPUT_DIR##*/}_combined_counts.csv ]] \
    && rm ${INPUT_DIR}/*.setmiss.annAB.counts.csv

exit $SUCCESS