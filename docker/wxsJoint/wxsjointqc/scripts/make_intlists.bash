#!/bin/bash
# This script sets up the intlists files.
INPUT_VCF=$1

NAMEBASE=${INPUT_VCF%.vcf.bgz}

# Make interval list from bed
$GATK BedToIntervalList \
	-I ${NAMEBASE}.bed \
	-O ${NAMEBASE}.interval_list \
	-SD ${REF_FASTA%.fasta}.dict

# Make directory to hold scattered intervals
mkdir ${INPUT_VCF%/*}/intlists

# Split to 50 intervals
$GATK SplitIntervals \
	-R ${REF_FASTA} \
	-L ${NAMEBASE}.interval_list \
	--scatter-count 50 \
	--interval-merging-rule OVERLAPPING_ONLY \
	--subdivision-mode INTERVAL_SUBDIVISION \
	-O ${INPUT_VCF%/*}/intlists

# Create beds from scattered interval lists for bcftools
for INT_LIST in ${INPUT_VCF%/*}/intlists/*.interval_list; do \
	$GATK IntervalListToBed \
		-I ${INT_LIST} \
		-O ${INT_LIST%.*}.bed
done
