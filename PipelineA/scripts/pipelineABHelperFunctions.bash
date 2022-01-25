# save argument to output directory.
function saveToOutputDirectory ()
{
	cp "$1"* ${OUTDIR}
}
function reportToLog ()
{
	echo -e "$1" >> ${LOGFILE}
}
# function for taking unmapped paired FASTQs to sorted BAM
function alignSortPairedFQs ()
{
bwa mem -M -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" ${REF_FASTA} $FQ1 $FQ2 | \
#${JAVA} -Djava.io.tmpdir=${WORKDIR} -Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false -jar ${PICARD} SortSam I=/dev/stdin O=${WORKDIR}/${SAMPLEID}.aln.srt.bam SO=coordinate CREATE_INDEX=true MAX_RECORDS_IN_RAM=2000000
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${OUTDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}/"
CURRENT_BAM="${OUTDIR}/${RGBASE}.aln.srt.bam"
samtools index -@ ${THREADS} -m "${MEM_SPLIT}G" ${CURRENT_BAM} ${CURRENT_BAM%.bam}
}
# function for taking unmapped interleaved FASTQs to sorted BAM
function alignSortInterleavedFQs ()
{
bwa mem -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" -M -p "${REF_FASTA}" "${INDIR}/${RGBASE}.fq.gz" | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${OUTDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}"
CURRENT_BAM="${OUTDIR}/${RGBASE}.aln.srt.bam"
samtools index -@ ${THREADS} -m "${MEM_SPLIT}G" ${CURRENT_BAM} ${CURRENT_BAM%.bam}
}
function getFreeMix ()
{
	/scripts/verifyBamID \
		--vcf "${REF_HAPMAP}" \
		--bam "${CURRENT_BAM}" \
		--chip-none \
		--maxDepth 1000 \
		--precise \
		--verbose \
		--ignoreRG \
		--out "${OUTDIR}/${RGBASE}_verifybam" \
		|& grep -v "Skipping marker"
	echo $(tail -n 1 ${OUTDIR}/${RGBASE}_verifybam.selfSM) | cut -f 6
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
bedtools intersect -u -a "$1" -b "$2" > ${OUTDIR}/${RGBASE}.isec.bam
CURRENT_BAM="${OUTDIR}/${RGBASE}.isec.bam"
}
# function for converting argument bam to cram
function saveBamAsCram ()
{
samtools view -T "${REF_FASTA}" -C -o "${OUTDIR}/${RGBASE}.cram" "$1"
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
java -Djava.io.tmpdir=${WORKDIR} -Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
-jar ${GATK360} -T DepthOfCoverage \
-R ${REF_FASTA} -nt 1 \
-ct 10 -ct 15 -ct 20 -ct 30 -ct 40 -ct 50 -ct 60 -ct 70 -ct 80 -ct 90 -ct 100 \
--omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable -L ${REF_PADBED} \
-I ${CURRENT_BAM} \
-o "${OUTDIR}/${SAMPLEID}_exome_coverage"
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
	-O "${OUTDIR}/${SAMPLEID}.recal.table1"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	ApplyBQSR \
	-R ${REF_FASTA} \
	-I ${CURRENT_BAM} \
	-bqsr-recal-file "${OUTDIR}/${SAMPLEID}.recal.table1" \
	-L ${REF_PADBED} \
	-O "${OUTDIR}/${SAMPLEID}.recal.bam"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	AnalyzeCovariates \
	-bqsr "${OUTDIR}/${SAMPLEID}.recal.table1" \
	-plots "${OUTDIR}/${SAMPLEID}_AnalyzeCovariates.pdf"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	BaseRecalibrator \
	-I ${CURRENT_BAM} \
	-R ${REF_FASTA} \
	--known-sites ${REF_MILLS_GOLD} \
	--known-sites ${REF_DBSNP} \
	--known-sites ${REF_ONEKGP1} \
	-O "${OUTDIR}/${SAMPLEID}.recal.table2"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	AnalyzeCovariates \
	-before "${OUTDIR}/${SAMPLEID}.recal.table1" \
	-after "${OUTDIR}/${SAMPLEID}.recal.table2" \
	-plots "${OUTDIR}/${SAMPLEID}_before-after-plots.pdf"
CURRENT_BAM="${OUTDIR}/${SAMPLEID}.recal.bam"
}
function callSampleVariants ()
{
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
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
	java -Djava.io.tmpdir=${WORKDIR} -Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
	-jar ${GATK360} -T VariantEval -R ${REF_FASTA} \
	-nt ${GATK_THREADS} \
	-L ${REF_PADBED} \
--dbsnp ${REF_DBSNP} \
--eval:"${SAMPLEID_VE}" "${CURRENT_VCF}" \
-o "${OUTDIR}/${SAMPLEID}_exome_varianteval.gatkreport"
}
