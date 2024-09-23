#!/bin/bash
CUR_VCF="$1"
MODE="$2"
RECAL_TABLE="$3"
NAMEBASE="${CUR_VCF%.*.*}.AS.${MODE}"
${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	ApplyVQSR \
	-AS \
	-R ${REF_FASTA} \
	-V ${CUR_VCF} \
	-mode ${MODE} \
	--truth-sensitivity-filter-level 99.7 \
	--recal-file ${RECAL_TABLE} \
	--tranches-file ${RECAL_TABLE%.*}.tranches \
	-O ${NAMEBASE}_recalibrated.vcf.gz \
&& echo "${NAMEBASE}_recalibrated.vcf.gz"
