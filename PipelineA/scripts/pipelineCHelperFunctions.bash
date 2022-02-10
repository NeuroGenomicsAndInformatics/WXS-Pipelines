function getCram ()
{
	rsync $STAGE_INDIR/$CRAM $INDIR/
}
function convertCramToBam ()
{
samtools view -b -@ $(($THREADS-1)) -T "${REF_FASTA}" -o "${INDIR}/${CRAM##*/}.test.bam" $CRAM
CURRENT_BAM="${INDIR}/${CRAM##*/}.test.bam"
samtools index -@ $(($THREADS-1)) $CURRENT_BAM
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
	-OUTPUT_DIR "${INDIR}" \
	-VALIDATION_STRINGENCY SILENT
}
function makeRgfiles ()
{
	SAMPLEID=$(echo $FULLSMID | cut -d^ -f1 | cut -d. -f1)
	PROJECTNAME=$(echo $FULLSMID | cut -d^ -f1 | cut -d.-f2 )
	BARCODE=$(echo $FULLSMID | cut -d^ -f2)
	PROJECT=$(echo $FULLSMID | cut -d^ -f3)
	local IFS=$'\n' RGS=($(samtools view -H "${CURRENT_BAM}" | grep "^@RG"))
	for RG in ${RGS[@]}; do
    RGID="$(echo ${RG} | grep -oP "(?<=ID:)[^[:space:]]*")"
    RGID_NEW="$(echo ${RGID} | cut -d: -f2- | sed 's/:/^/g')"
    if [ -f "${INDIR}/${RGID//:/_}_2.fastq" ]; then
		mv -vf "${INDIR}/${RGID//:/_}_1.fastq" "${INDIR}/${FULLSMID}.${RGID_NEW}.r1.fastq"
		mv -vf "${INDIR}/${RGID//:/_}_2.fastq" "${INDIR}/${FULLSMID}.${RGID_NEW}.r2.fastq";
	 	else
		mv -vf "${INDIR}/${RGID//:/_}_1.fastq" "${INDIR}/${FULLSMID}.${RGID_NEW}.fastq"
		fi
		echo "@RG\tID:${RGID}\tPL:illumina\tPU:${RGID}:${BARCODE}\tLB:${BARCODE}\tSM:${SAMPLEID}.${PROJECTNAME}\tDS:${FULLSMID}" > "${INDIR}/${FULLSMID}.${RGID_NEW}.rgfile"
  done
}
