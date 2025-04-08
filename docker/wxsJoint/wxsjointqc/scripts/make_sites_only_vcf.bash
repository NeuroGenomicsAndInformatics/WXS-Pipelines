#!/bin/bash
# This script makes a sites-only vcf.
INPUT_DIR=$1
INTERVAL=$2
[[ -z $INTERVAL ]] && INTERVAL=${LSB_JOBINDEX}
INPUT_VCF=$(find ${INPUT_DIR} -name "*.${INTERVAL}.*.splitMA.geno.vcf.gz")
OUTPUT_VCF="${INPUT_VCF%.*.*}.sites.vcf.gz"

$GATK MakeSitesOnlyVcf \
    -I $INPUT_VCF \
    -O $OUTPUT_VCF
    