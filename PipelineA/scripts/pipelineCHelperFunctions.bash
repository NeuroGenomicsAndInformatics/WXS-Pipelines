# save argument to output directory.
function saveToOutputDirectory ()
{
	cp "$1*" ${OUTDIR}
}
function reportToLog ()
{
	echo -e "$1" >> ${LOGFILE}
}
# function for converting argument bam to cram
function saveBamAsCram ()
{
samtools view -T "${REF_FASTA}" -C -o "${OUTDIR}/${RGBASE}.cram" "$1"
}
function convertCramToBam ()
{
samtools view -b -T "${REF_FASTA}" -o "${WORKDIR}/${RGBASE}.bam" "$1"
CURRENT_BAM="${WORKDIR}/${RGBASE}.bam"
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
function extractRGs ()
{
	IFS=$'\n' RGS=($(samtools view -H "${CURRENT_BAM}" | grep "^@RG"))
	for RG in ${RGS[@]}; do RGID="$(echo ${RG} | grep -oP "(?<=ID:)[^[:space:]]*")"
	RGID_NEW="$(echo ${RGID} | cut -d: -f2- | sed 's/:/^/g')"
	RGBASE="${FULLSM}.${RGID_NEW}"
}
