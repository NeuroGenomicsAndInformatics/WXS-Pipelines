#!/bin/bash
JOINT_VCF="$1"
REGION="$2"
NAMEBASE="${JOINT_VCF%.*.*}"
bcftools norm \
	--threads $LSB_MAX_NUM_PROCESSORS \
	-r $REGION \
	-m -any \
	-o ${NAMEBASE}-norm.vcf.gz \
	${JOINT_VCF} \
	&& tabix -s1 -b2 -e2 ${NAMEBASE}-norm.vcf.gz
bcftools norm \
	--threads $LSB_MAX_NUM_PROCESSORS \
	-r $REGION \
	-f /scratch1/fs1/cruchagac/WXSref/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
	-o ${NAMEBASE}-norm-ref.vcf.gz \
	${NAMEBASE}-norm.vcf.gz \
	&& tabix -s1 -b2 -e2 ${NAMEBASE}-norm-ref.vcf.gz
bcftools annotate \
	--threads $LSB_MAX_NUM_PROCESSORS \
	-r $REGION \
	-x ID -I +'%CHROM:%POS:%REF:%ALT' \
	-o ${NAMEBASE}-split.vcf.gz \
	${NAMEBASE}-norm-ref.vcf.gz \
	&& tabix -s1 -b2 -e2 ${NAMEBASE}-split.vcf.gz
rm ${NAMEBASE}-norm.vcf.gz
rm ${NAMEBASE}-norm-ref.vcf.gz