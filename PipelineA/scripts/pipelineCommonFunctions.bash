function stageDataForRGBASE ()
{
	rsync -rg  $STAGE_DIR/${RGBASE}* $INDIR
}
# save argument to output directory.
function saveToOutputDirectory ()
{
	rsync "$1"* ${OUTDIR}
	if [[ -e "${1::-1}i"]]; then cp "${1::-1}i" ${OUTDIR}; fi
}
# save argument to final output directory.
function saveToFinalOutputDirectory ()
{
	rsync "$1"* ${FINAL_OUTDIR}
	if [[ -e "${1::-1}i"]; then cp "${1::-1}i" ${FINAL_OUTDIR}; fi
}
# prints argument to stdin. could consider for stderr
function reportToLog ()
{
	echo -e "\n\n ===================================== \n\n"
	echo -e "$1"
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
# function for getting the size of the file given as an argument
function getFileSize ()
{
	echo $(wc -c "$1" | cut -d' ' -f1)
}
