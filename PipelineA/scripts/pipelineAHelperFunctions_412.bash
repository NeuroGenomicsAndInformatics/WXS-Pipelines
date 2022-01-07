# save argument to output directory.
function saveToOutputDirectory ()
{
	cp "$1" ${OUTDIR}
}
function chmodWorkDirectory ()
{
	chmod -r 777 ${WOrKDIR}
}
function reportToLog ()
{
	echo "$1" >> ${LOGFILE}
}
# function for taking unmapped FASTQs to sorted BAM
function alignSortFQsAndVerifyBam ()
{
bwa mem -M -t ${THREADS} -R "${RG}" ${REF_FASTA} ${FQ1} ${FQ2} | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o "${WORKDIR}/${FULLSMID}.aln.srt.bam" -T "${WORKDIR}/"
CURRENT_BAM="${WORKDIR}/${FULLSMID}.aln.srt.bam"
#/scripts/verifyBamID \
#	--vcf "${REF_HAPMAP}" \
#	--bam "${CURRENT_BAM}" \
#	--chip-none \
#	--maxDepth 1000 \
#	--precise \
#	--verbose \
#	--ignoreRG \
#	--out "${OUTDIR}/${FULLSMID}_verifybam" \
#	|& grep -v "Skipping marker"
chmodWorkDirectory
}
function getFreeMix ()
{
	echo $(tail -n 1 ${OUTDIR}/${FULLSMID}_verifybam.selfSM) | cut -f 6
}
function validateCurrentBam ()
{
java -Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
-jar ${PICARD} \
ValidateSamFile -I ${CURRENT_BAM} \
-R ${REF_FASTA} \
--CREATE_INDEX true \
--IGNORE MISSING_READ_GROUP \
--IGNORE RECORD_MISSING_READ_GROUP \
--IGNORE INVALID_VERSION_NUMBER \
--IGNORE INVALID_TAG_NM \
--TMP_DIR ${WORKDIR} \
-MODE SUMMARY
}
# first argument is for the bam file and second is for the reference bed
function intersectBamWithBed ()
{
bedtools intersect -u -a "$1" -b "$2" > "${WORKDIR}/${FULLSMID}.isec.bam"
CURRENT_BAM="${WORKDIR}/${FULLSMID}.isec.bam"
}
# function for converting argument bam to cram
function saveBamAsCram ()
{
samtools view -T "${REF_FASTA}" -C -o "${OUTDIR}/${FULLSMID}.cram" "$1"
}
function markDuplicates ()
{
java -Djava.io.tmpdir=${WORKDIR} \
-jar ${PICARD} MarkDuplicates -I ${CURRENT_BAM} \
-O "${WORKDIR}/${FULLSMID}_GATKready.bam" \
-M "${WORKDIR}/${FULLSMID}_metrics.txt" \
-QUIET true \
-MAX_RECORDS_IN_RAM 2000000 \
-ASSUME_SORTED TRUE \
-CREATE_INDEX TRUE
CURRENT_BAM="${WORKDIR}/${FULLSMID}_GATKready.bam"
}
function analyzeDepthOfCoverage ()
{
${GATK} --java-options "-Djava.io.tmpdir=${WORKDIR} -Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false"\
DepthOfCoverage \
-R ${REF_FASTA} -nt 1 \
-ct 10 -ct 15 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 \
--omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable -L ${REF_COVBED} \
-I ${CURRENT_BAM} \
-o "${WORKDIR}/${FULLSMID}_exome_coverage"
}
# function for base recalibration and plots
function recalibrateBases ()
{
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	BaseRecalibrator \
	-I ${CURRENT_BAM} \
	-R ${REF_FASTA} \
	--known-sites ${REF_MILLS_GOLD} \
	--known-sites ${REF_DBSNP} \
	--known-sites ${REF_ONEKGP1} \
	-O "${WORKDIR}/${FULLSMID}.recal.table1"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	ApplyBQSR \
	-R ${REF_FASTA} \
	-I ${CURRENT_BAM} \
	-bqsr-recal-file "${WORKDIR}/${SAMPLEID}.recal.table1" \
	-L ${REF_COVBED} \
	-O "${WORKDIR}/${FULLSMID}.recal.bam"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	AnalyzeCovariates \
	-bqsr "${WORKDIR}/${FULLSMID}.recal.table1" \
	-plots "${WORKDIR}/${FULLSMID}_AnalyzeCovariates.pdf"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	BaseRecalibrator \
	-I ${CURRENT_BAM} \
	-R ${REF_FASTA} \
	--known-sites ${REF_MILLS_GOLD} \
	--known-sites ${REF_DBSNP} \
	--known-sites ${REF_ONEKGP1} \
	-O "${WORKDIR}/${FULLSMID}.recal.table2"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	AnalyzeCovariates \
	-before "${WORKDIR}/${FULLSMID}.recal.table1" \
	-after "${WORKDIR}/${FULLSMID}.recal.table2" \
	-plots "${WORKDIR}/${FULLSMID}_before-after-plots.pdf"
CURRENT_BAM="${WORKDIR}/${FULLSMID}.recal.bam"
}
function callSampleVariants ()
{
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
HaplotypeCaller -R ${REF_FASTA} \
	-I ${CURRENT_BAM} \
	-ERC GVCF \
	--dbsnp ${REF_DBSNP} \
	-L ${REF_COVBED} \
	-O "${WORKDIR}/${FULLSMID}.raw.snps.indels.g.vcf.gz"
CURRENT_VCF="${WORKDIR}/${FULLSMID}.raw.snps.indels.g.vcf.gz"
}
function evaluateSampleVariants ()
{
${GATK} VariantEval -R ${REF_FASTA} \
-L ${REF_CCDS} -nt ${GATK_THREADS} \
--dbsnp ${REF_DBSNP} \
--eval ${SAMPLEID_VE} "${WORKDIR}/${FULLSMID}.raw.snps.indels.g.vcf.gz" -o "${WORKDIR}/${FULLSMID}_exome_varianteval.gatkreport"
}
