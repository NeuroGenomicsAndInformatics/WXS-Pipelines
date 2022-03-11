function stageDataForSample ()
{
	rsync ${STAGE_INDIR}/*.ba* ${INDIR}
	rsync ${STAGE_INDIR}/stage1complete.txt ${INDIR}
}
function transferOutputFilesToStorage ()
{
	mkdir ${FINAL_OUTDIR}
	rsync ${OUTDIR}/*.env ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*.rgfile ${FINAL_OUTDIR}/
	rsync ${STAGE_INDIR}/*.cram ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*_GATK* ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*vcf* ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*verifybam* ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*_exome_* ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*.pdf ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*.recal.table* ${FINAL_OUTDIR}/
}
function cleanUp ()
{
	rsync ${OUTDIR}/*.out /output/logs/${FULLSMID}
	rsync ${OUTDIR}/*.txt /output/logs/${FULLSMID}
	rm -R ${OUTDIR}
	rm -R ${INDIR}
	rm -d ${STAGE_INDIR}
}
function markDuplicates ()
{
java -Djava.io.tmpdir=${WORKDIR} \
-jar ${PICARD} MarkDuplicates ${MD_INPUTS[@]} \
O="${OUTDIR}/${SAMPLEID}_GATKready.bam" \
M="${OUTDIR}/${SAMPLEID}_metrics.txt" \
QUIET=true \
MAX_RECORDS_IN_RAM=2000000 \
ASSUME_SORTED=TRUE \
CREATE_INDEX=TRUE \
TMP_DIR=${WORKDIR}
CURRENT_BAM="${OUTDIR}/${SAMPLEID}_GATKready.bam"
}
function analyzeDepthOfCoverage ()
{
java -Djava.io.tmpdir=${WORKDIR} -Xms2g -Xmx${S2MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
-jar ${GATK360} -T DepthOfCoverage \
-R ${REF_FASTA} -nt 1 \
-ct 10 -ct 15 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 \
--omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable -L ${REF_PADBED} \
-I ${CURRENT_BAM} \
-o "${OUTDIR}/${SAMPLEID}_exome_coverage"
}
function getAverageDepthOfCoverage ()
{
	tail -n1 ${OUTDIR}/${SAMPLEID}_exome_coverage.sample_summary | cut -f3
}
# function for base recalibration and plots
function recalibrateBases ()
{
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S2MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	BaseRecalibrator \
	-I ${CURRENT_BAM} \
	-R ${REF_FASTA} \
	--known-sites ${REF_MILLS_GOLD} \
	--known-sites ${REF_DBSNP} \
	--known-sites ${REF_ONEKGP1} \
	-O "${OUTDIR}/${SAMPLEID}.recal.table1"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S2MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	ApplyBQSR \
	-R ${REF_FASTA} \
	-I ${CURRENT_BAM} \
	-bqsr-recal-file "${OUTDIR}/${SAMPLEID}.recal.table1" \
	-L ${REF_PADBED} \
	-O "${OUTDIR}/${SAMPLEID}.recal.bam"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S2MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	AnalyzeCovariates \
	-bqsr "${OUTDIR}/${SAMPLEID}.recal.table1" \
	-plots "${OUTDIR}/${SAMPLEID}_AnalyzeCovariates.pdf"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S2MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	BaseRecalibrator \
	-I ${CURRENT_BAM} \
	-R ${REF_FASTA} \
	--known-sites ${REF_MILLS_GOLD} \
	--known-sites ${REF_DBSNP} \
	--known-sites ${REF_ONEKGP1} \
	-O "${OUTDIR}/${SAMPLEID}.recal.table2"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S2MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	AnalyzeCovariates \
	-before "${OUTDIR}/${SAMPLEID}.recal.table1" \
	-after "${OUTDIR}/${SAMPLEID}.recal.table2" \
	-plots "${OUTDIR}/${SAMPLEID}_before-after-plots.pdf"
CURRENT_BAM="${OUTDIR}/${SAMPLEID}.recal.bam"
}
function verifyBamID ()
{
	/scripts/verifyBamID \
		--vcf "${REF_HAPMAP}" \
		--bam "${CURRENT_BAM}" \
		--chip-none \
		--maxDepth 1000 \
		--precise \
		--verbose \
		--ignoreRG \
		--out "${OUTDIR}/${FULLSMID}_verifybam" \
		|& grep -v "Skipping marker"
}
function getFreeMix ()
{
	tail -n1 ${OUTDIR}/${FULLSMID}_verifybam.selfSM | cut -f7
}

function callSampleVariants ()
{
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S2MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
HaplotypeCaller -R ${REF_FASTA} \
	-I ${CURRENT_BAM} \
	-ERC GVCF \
	--dbsnp ${REF_DBSNP} \
	-L ${REF_PADBED} \
	-O "${OUTDIR}/${SAMPLEID}.raw.snps.indels.g.vcf.gz"
CURRENT_VCF="${OUTDIR}/${SAMPLEID}.raw.snps.indels.g.vcf.gz"
}
function evaluateSampleVariants ()
{
	java -Djava.io.tmpdir=${WORKDIR} -Xms2g -Xmx${S2MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
	-jar ${GATK360} -T VariantEval -R ${REF_FASTA} \
	-nt ${GATK_THREADS} \
	-L ${REF_PADBED} \
--dbsnp ${REF_DBSNP} \
--eval:"${SAMPLEID_VE}" "${CURRENT_VCF}" \
-o "${OUTDIR}/${SAMPLEID}_exome_varianteval.gatkreport"
}
function getTitvRatio ()
{
	echo $(cat ${OUTDIR}/${SAMPLEID}_exome_varianteval.gatkreport | grep TiTvVariantEvaluator | grep all) | cut -d ' ' -f 8
}
