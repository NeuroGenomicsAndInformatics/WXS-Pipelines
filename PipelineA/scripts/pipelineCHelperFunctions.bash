# save argument to output directory.
function saveToOutputDirectory ()
{
	cp "$1" ${OUTDIR}
}
function reportToLog ()
{
	echo -e "$1" >> ${LOGFILE}
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
# function for converting argument bam to cram
function saveBamAsCram ()
{
samtools view -T "${REF_FASTA}" -C -o "${OUTDIR}/${RGBASE}.cram" "$1"
}
function revertBamToFastqs ()
{
	java -Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
	-jar ${PICARD} \
	RevertSam -I ${CURRENT_BAM} \
	-O /dev/stdout \
	-SORT_ORDER queryname \
	-COMPRESSION_LEVEL 0 \
	-VALIDATION_STRINGENCY SILENT \
	| java -Xms2g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
	-jar ${PICARD} \
	SamToFastq -I /dev/stdin \
	-OUTPUT_PER_RG true \
	-RG_TAG ID \
	-OUTPUT_DIR "${OUTDIR}" \
	-VALIDATION_STRINGENCY SILENT
}
