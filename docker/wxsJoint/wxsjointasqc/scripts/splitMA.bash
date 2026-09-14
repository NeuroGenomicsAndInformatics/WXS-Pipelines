#!/bin/bash
JOINT_VCF="$1"
NAMEBASE="${JOINT_VCF%.*.*}"
bcftools norm \
	--threads $LSB_MAX_NUM_PROCESSORS \
	-m -any \
	-o ${NAMEBASE}.splitMA.vcf.gz \
	${JOINT_VCF} \
	&& tabix -s1 -b2 -e2 ${NAMEBASE}.splitMA.vcf.gz
