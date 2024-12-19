#!/bin/bash
MODE="INDEL"
NAMEBASE="${JOINT_VCF%.*.*}.AS.${CHR}.${MODE}"
${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	VariantRecalibrator \
	-AS \
	-R ${REF_FASTA} \
	-L ${CHR} \
	-V ${JOINT_VCF} \
	--resource:mills,known=false,training=true,truth=true,prior=12.0 ${REF_MILLS_GOLD} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${REF_DBSNP} \
	-an AS_QD -an AS_DP -an AS_MQRankSum -an AS_ReadPosRankSum -an AS_FS -an AS_SOR \
	-mode $MODE \
	--max-gaussians 4 \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.0 -tranche 90.0 \
	-O ${NAMEBASE}_recalibrate.recal \
	--tranches-file ${NAMEBASE}_recalibrate.tranches \
	--rscript-file ${NAMEBASE}_recalibrate_plots.R \
&& echo "${NAMEBASE}_recalibrate.recal"
