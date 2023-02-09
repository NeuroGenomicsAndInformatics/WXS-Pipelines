#!/bin/bash
MODE="SNP"
NAMEBASE="${JOINT_VCF%.*.*}.AS.${CHR}.${MODE}"
/ref/gatk-4.2.6.1_mod/gatk \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	VariantRecalibrator \
	-AS \
	-R ${REF_FASTA} \
	-L ${CHR} \
	-V ${JOINT_VCF} \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${REF_HAPMAP} \
	--resource:omni,known=false,training=true,truth=true,prior=12.0 /ref/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${REF_ONEKGP1} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${REF_DBSNP} \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
	-mode $MODE \
	--max-gaussians 2 \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.0 -tranche 90.0 \
	-O ${NAMEBASE}_recalibrate.recal \
	--tranches-file ${NAMEBASE}_recalibrate.tranches \
	--rscript-file ${NAMEBASE}_recalibrate_plots.R \
&& echo "${NAMEBASE}_recalibrate.recal"
