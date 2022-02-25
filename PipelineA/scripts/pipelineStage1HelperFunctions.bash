function stageDataForRGBASE ()
{
	rsync -L ${STAGE_INDIR}/${RGBASE}* ${INDIR}
	rsync ${INDIR}/${RGBASE}.rgfile ${OUTDIR}
}
function cleanUp ()
{
	rm ${INDIR}/${RGBASE}*
}
# function for taking unmapped paired FASTQs to sorted BAM
function alignSortPairedFQs ()
{
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
bwa mem -M -t ${S1THREADS} -R $(head -n1 ${INDIR}/${RGBASE}.rgfile) ${REF_FASTA} $FQ1 $FQ2 | \
samtools sort -@ $((${S1THREADS} / 2)) -m "$((${MEM_SPLIT} * 2))G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}/"
}
function alignSortPairedHugeFQs ()
{
CURRENT_BAM="${OUTDIR}/${RGBASE}.aln.srt.bam"
bwa mem -M -t ${S1THREADS} -R $(head -n1 ${INDIR}/${RGBASE}.rgfile) ${REF_FASTA} $FQ1 $FQ2 | \
samtools view -b -1 -o ${OUTDIR}/${RGBASE}.aln.bam \
&& samtools sort -@ $((${S1THREADS} / 2)) -m "$((${MEM_SPLIT} * 2))G" -o ${OUTDIR}/${RGBASE}.aln.srt.bam -T "${OUTDIR}/" "${OUTDIR}/${RGBASE}.aln.bam"
}
# function for taking unmapped interleaved FASTQs to sorted BAM
function alignSortInterleavedFQs ()
{
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
bwa mem -t ${S1THREADS} -R $(head -n1 ${INDIR}/${RGBASE}.rgfile) -M -p "${REF_FASTA}" $FQI | \
samtools sort -@ ${S1THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}"
}
function alignSortHugeInterleavedFQs ()
{
CURRENT_BAM="${WORKDIR}/${RGBASE}.aln.srt.bam"
bwa mem -t ${S1THREADS} -R $(head -n1 ${INDIR}/${RGBASE}.rgfile) -M -p "${REF_FASTA}" "${INDIR}/${RGBASE}.fq.gz" | \
samtools view -b -1 -o ${WORKDIR}/${RGBASE}.aln.bam \
&& samtools sort -@ ${S1THREADS} -m "${MEM_SPLIT}G" -o ${WORKDIR}/${RGBASE}.aln.srt.bam -T "${WORKDIR}/" "${WORKDIR}/${RGBASE}.aln.bam"
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
