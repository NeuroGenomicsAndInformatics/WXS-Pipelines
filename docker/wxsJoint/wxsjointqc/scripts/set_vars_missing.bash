#!/bin/bash
# This script sets genotypes to missing in a vcf file.
INPUT_VCF=$1
INTERVAL=$2
[[ -z $INTERVAL ]] && INTERVAL=${LSB_JOBINDEX}
INT_BEDS=($(ls ${INPUT_VCF%/*}/intlists | grep scattered.bed))
INT_BED="${INPUT_VCF%/*}/intlists/${INT_BEDS[${INTERVAL}-1]}"
OUTPUT_VCF="${INPUT_VCF%.*.*}.${INTERVAL}.setmiss.vcf.gz"

bcftools +setGT \
	-R ${INT_BED} \
	-Ou \
	${INPUT_VCF} \
	-- \
	-t q \
	-n . \
	-i 'FORMAT/DP<10' \
| bcftools +setGT \
	-Ou \
	-- \
	-t q \
	-n . \
	-i 'FORMAT/GQ<20' \
| bcftools annotate \
	--set-id '%CHROM\:%POS\:%REF\:%ALT' \
	-Oz \
	-o $OUTPUT_VCF \
	-W=tbi
