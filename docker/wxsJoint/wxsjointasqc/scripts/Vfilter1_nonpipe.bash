#!/bin/bash
RECAL_VCF="$1"
NAMEBASE="${RECAL_VCF%.*.*}"
COUNT_FILE="${NAMEBASE}.counts.csv"

echo "Filter,File Location,Number Variants" > $COUNT_FILE

echo "Recal,${NAMEBASE}.vcf.gz,$(zcat ${NAMEBASE}.vcf.gz | grep -v '^#' | wc -l)" > $COUNT_FILE
echo "${NAMEBASE}.vcf.gz" > ${NAMEBASE}.vcftools.log; zcat ${NAMEBASE}.vcf.gz | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${NAMEBASE}.vcftools.log

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${INT_LIST} \
	-V ${NAMEBASE}.vcf.gz \
	--exclude-filtered \
	-O ${NAMEBASE}-PASS.vcf.gz \
&& rm ${RECAL_VCF} && ${RECAL_VCF/INDEL/SNP}

echo "Filtered,${NAMEBASE}-PASS.vcf.gz,$(zcat ${NAMEBASE}-PASS.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE
echo "${NAMEBASE}-PASS.vcf.gz" >> ${NAMEBASE}-PASS.vcftools.log; zcat ${NAMEBASE}-PASS.vcf.gz | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${NAMEBASE}-PASS.vcftools.log

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	VariantsToTable \
	-R ${REF_FASTA} \
	-L ${INT_LIST} \
	-V ${NAMEBASE}-PASS.vcf.gz \
	--show-filtered \
	--split-multi-allelic \
	-F CHROM -F POS -F DP \
	-O ${NAMEBASE}-PASS-DP.table

/scripts/CALCULATE_DP.R ${NAMEBASE}-PASS-DP.table

SNV_DP_TR=`tail -n +2 ${NAMEBASE}-PASS-DP.calc | cut -f 3`

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${INT_LIST} \
	-V ${NAMEBASE}-PASS.vcf.gz \
	-select "DP <= ${SNV_DP_TR}" \
	-O ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}.vcf.gz \
&& rm ${NAMEBASE}-PASS.vcf.gz

echo "DP Filtered,${NAMEBASE}-PASS-maxDP${SNV_DP_TR}.vcf.gz,$(zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE
echo "${NAMEBASE}-PASS-maxDP${SNV_DP_TR}.vcf.gz" >> ${NAMEBASE}-PASS.DP.vcftools.log; zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}.vcf.gz | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${NAMEBASE}-PASS.DP.vcftools.log
	
LCR="/scripts/LCR-hs38.bed"
bedtools subtract -header -A -a ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}.vcf.gz -b ${LCR} | bgzip -c > ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz 
bcftools index -t --threads $LSB_MAX_NUM_PROCESSORS ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz \
&& rm ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}.vcf.gz

echo "LCR Filtered,${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz,$(zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz | grep -v '#' | wc -l)" >> $COUNT_FILE
echo "${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz" >> ${NAMEBASE}-PASS.DP.LCR.vcftools.log; zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${NAMEBASE}-PASS.DP.LCR.vcftools.log

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${INT_LIST} \
	-V ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz \
	--exclude-non-variants \
	-O ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants.vcf.gz \
&& rm ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR.vcf.gz

echo "Non-Variant Filtered,${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants.vcf.gz,$(zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE
echo "${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants.vcf.gz" >> ${NAMEBASE}-PASS.DP.LCR.nonVar.vcftools.log; zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants.vcf.gz | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${NAMEBASE}-PASS.DP.LCR.nonVar.vcftools.log

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${INT_LIST} \
	-V ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants.vcf.gz \
	-select 'vc.hasGenotypes() && vc.getCalledChrCount(vc.getAltAlleleWithHighestAlleleCount())/(1.0*vc.getCalledChrCount()) < 1.0'  \
	-O ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1.vcf.gz \
&& rm ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants.vcf.gz

echo "AF Filtered,${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1.vcf.gz,$(zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE
echo "${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1.vcf.gz" >> ${NAMEBASE}-PASS.DP.LCR.nonVar.AF1.vcftools.log; zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1.vcf.gz | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${NAMEBASE}-PASS.DP.LCR.nonVar.AF1.vcftools.log

java -Xmx80g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -XX:+UseSerialGC \
	-jar /ref/GATK360.jar \
	-T VariantAnnotator \
	-A AlleleBalance \
	-R /ref/20190812_GATK_38_googlebundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
	-V ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1.vcf.gz \
	-o ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABannotated.vcf.gz \
&& rm ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1.vcf.gz

${GATK} \
    --java-options "-Xmx80g -XX:ConcGCThreads=1 -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	SelectVariants \
	-R ${REF_FASTA} \
	-L ${INT_LIST} \
	-V ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABannotated.vcf.gz \
	-select '(!vc.hasAttribute("ABHet")) || (ABHet >= 0.30 && ABHet <= 0.70)' \
	-O ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABfiltered.vcf.gz
SUCCESS=$?

[[ $SUCCESS -eq 0 ]] && rm ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABannotated.vcf.gz

echo "AB Filtered,${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABfiltered.vcf.gz,$(zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABfiltered.vcf.gz | grep -v '^#' | wc -l)" >> $COUNT_FILE
echo "${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABfiltered.vcf.gz" >> ${NAMEBASE}-PASS.DP.LCR.nonVar.AF1.ABHet.vcftools.log; zcat ${NAMEBASE}-PASS-maxDP${SNV_DP_TR}-LCR-nonVariants-AF1-ABfiltered.vcf.gz | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${NAMEBASE}-PASS.DP.LCR.nonVar.AF1.ABHet.vcftools.log

exit $SUCCESS