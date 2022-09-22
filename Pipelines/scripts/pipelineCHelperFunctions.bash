function getCram ()
{
	rsync -L $STAGE_INDIR/$CRAM $INDIR/
}
function convertCramToBam ()
{
CURRENT_BAM="${INDIR}/${CRAM}.bam"
samtools view -b -@ $(($THREADS-1)) -T "${REF_FASTA}" -o "${CURRENT_BAM}" ${INDIR}/${CRAM}
samtools index -@ $(($THREADS-1)) $CURRENT_BAM
}
function revertBamToFastqs ()
{
	java -Xms${MEM_SPLIT}g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
	-jar ${PICARD} \
	RevertSam -I ${CURRENT_BAM} \
	-O /dev/stdout \
	-SORT_ORDER queryname \
	-COMPRESSION_LEVEL 0 \
	-VALIDATION_STRINGENCY SILENT \
	--TMP_DIR ${INDIR} \
	| java -Xms${MEM_SPLIT}g -Xmx${MEM}g -XX:+UseSerialGC -Dpicard.useLegacyParser=false \
	-jar ${PICARD} \
	SamToFastq -I /dev/stdin \
	-OUTPUT_PER_RG true \
	-RG_TAG ID \
	-OUTPUT_DIR "${INDIR}" \
	-VALIDATION_STRINGENCY SILENT \
	--TMP_DIR ${INDIR}
}
function makeRgfiles ()
{
	SAMPLEID=$(echo $FULLSMID | cut -d^ -f1)
	BARCODE=$(echo $FULLSMID | cut -d^ -f2)
	PROJECT=$(echo $FULLSMID | cut -d^ -f3)
	samtools view -H ${CURRENT_BAM} | grep "^@RG" > ${WORKDIR}/RGS.txt
	cat ${WORKDIR}/RGS.txt | while read RG; do
    RGID="$(echo ${RG} | grep -oP "(?<=ID:)[^[:space:]]*")"
		RGID_NEW="$(echo ${RG} | grep -oP "(?<=PU:)[^[:space:]]*" | sed 's/./^/g')"
    if [ -f "${INDIR}/${RGID}_2.fastq" ]; then
		mv -vf "${INDIR}/${RGID}_1.fastq" "${INDIR}/${FULLSMID}.${RGID_NEW}.r1.fastq"
		mv -vf "${INDIR}/${RGID}_2.fastq" "${INDIR}/${FULLSMID}.${RGID_NEW}.r2.fastq";
	 	else
		mv -vf "${INDIR}/${RGID}_1.fastq" "${INDIR}/${FULLSMID}.${RGID_NEW}.fastq"
		fi
		echo "@RG\tID:${RGID_NEW//^/:}\tPL:illumina\tPU:${RGID_NEW//^/:}:${BARCODE}\tLB:${BARCODE}\tSM:${SAMPLEID}\tDS:${FULLSMID}" > ${INDIR}/${FULLSMID}.${RGID_NEW}.rgfile
  done
}
function cleanUp ()
{
	rsync ${INDIR}/*.rgfile ${STAGE_INDIR}/
	rsync ${INDIR}/*.fastq ${STAGE_INDIR}/
	rm -R ${INDIR}
}
