function stageDataForRGBASE ()
{
	rsync -rg  $STAGE_DIR/${RGBASE}* $INDIR
}
# function for taking unmapped paired FASTQs to sorted BAM
function alignSortPairedFQs ()
{
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
bwa mem -M -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" ${REF_FASTA} $FQ1 $FQ2 | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}/")
}
function alignSortPairedHugeFQs ()
{
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
bwa mem -M -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" ${REF_FASTA} $FQ1 $FQ2 | \
samtools view -b -1 -o ${WORKDIR}/${RGBASE}.aln.bam \
&& samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}/" "${WORKDIR}/${RGBASE}.aln.bam"
}
# function for taking unmapped interleaved FASTQs to sorted BAM
function alignSortInterleavedFQs ()
{
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
bwa mem -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" -M -p "${REF_FASTA}" $FQI | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}")
}
function alignSortHugeInterleavedFQs ()
{
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
bwa mem -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" -M -p "${REF_FASTA}" "${INDIR}/${RGBASE}.fq.gz" | \
samtools view -b -1 -o ${WORKDIR}/${RGBASE}.aln.bam \
&& samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}/" "${WORKDIR}/${RGBASE}.aln.bam"
}
# first argument is for the bam file and second is for the reference bed
function intersectBamWithBed ()
{
bedtools intersect -u -a "$1" -b "$2" > ${WORKDIR}/${RGBASE}.isec.bam
CURRENT_BAM="${WORKDIR}/${RGBASE}.isec.bam"
}
# function for converting argument bam to cram
function saveBamAsCram ()
{
samtools view -T "${REF_FASTA}" -C -o "${OUTDIR}/${RGBASE}.cram" "$1"
}
