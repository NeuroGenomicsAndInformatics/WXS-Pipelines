function stageDataForRGBASE ()
{
	rsync -rg  $STAGE_DIR/${RGBASE}* $INDIR
}
# function for taking unmapped paired FASTQs to sorted BAM
function alignSortPairedFQs ()
{
bwa mem -M -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" ${REF_FASTA} $FQ1 $FQ2 | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}/"
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
}
function alignSortPairedHugeFQs ()
{
bwa mem -M -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" ${REF_FASTA} $FQ1 $FQ2 | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${OUTDIR}/${RGBASE}.aln.srt.bam -T "${OUTDIR}/"
CURRENT_BAM="${OUTDIR}/${RGBASE}.aln.srt.bam"
}
# function for taking unmapped interleaved FASTQs to sorted BAM
function alignSortInterleavedFQs ()
{
bwa mem -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" -M -p "${REF_FASTA}" $FQI | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}"
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
}
function alignSortHugeInterleavedFQs ()
{
bwa mem -t ${THREADS} -R "$(<${INDIR}/${RGBASE}.rgfile)" -M -p "${REF_FASTA}" "${INDIR}/${RGBASE}.fq.gz" | \
samtools sort -@ ${THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${OUTDIR}"
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
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
