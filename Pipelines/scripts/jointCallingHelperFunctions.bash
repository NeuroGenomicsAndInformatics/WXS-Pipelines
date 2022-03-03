function stageDataForCOHORT ()
{
	rsync -rpL $STAGE_INDIR/ $INDIR
}
function transferOutputFilesToStorage ()
{
	mkdir ${FINAL_OUTDIR}
	rsync ${OUTDIR}/*.env ${FINAL_OUTDIR}/
	rsync ${OUTDIR}/*vcf* ${FINAL_OUTDIR}/
	#rm -R ${OUTDIR}
	#rm -R ${INDIR}
}
function cleanUp ()
{
	#rm -R $STAGE_INDIR
	rm -R $INDIR
	#rm -R $OUTDIR
}
function makeSampleMap ()
{
	if [[ -r ${INDIR}/${SAMPLE_MAP##*/} ]]; then
	SAMPLE_MAP=${${INDIR}/${SAMPLE_MAP##*/}}
	else
	SAMPLE_MAP="${WORKDIR}/SampleMap.txt"
	cp ${SAMPLE_MAP} ${WORKDIR}/SampleMap.txt || \
 	(echo -n "" > ${SAMPLE_MAP}
	for SAMPLE in $(ls -d ${INDIR}) ; do
	echo -en "${echo $SAMPLE | cut -d'^' -f1}\t" >> ${SAMPLE_MAP}
	echo "$(find ${STAGE_DIR}/${SAMPLE} -name "*.vcf.gz" -print)" >> ${SAMPLE_MAP}
done)
fi
}
function buildGenomicDB ()
{
	DATABASE="${WORKDIR}/db"
if [[ ! -d ${OUTDIR}/dbwork ]]; then mkdir ${OUTDIR}/dbwork; fi
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S3MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	GenomicsDBImport \
	${INTERVAL} \
	--sample-name-map ${SAMPLE_MAP} \
	--genomicsdb-workspace-path ${WORKDIR}/db \
	--batch-size 50 \
	--tmp-dir $OUTDIR/dbwork \
	--reader-threads ${S3THREADS}
}
function jointCallCohort ()
{
	JOINT_GVCF="${OUTDIR}/${COHORT}.joint.g.vcf.gz"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S3MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
GenotypeGVCFs \
	-R ${REF_FASTA} \
	-V gendb://${DATABASE} \
	-O "${OUTDIR}/${COHORT}.joint.g.vcf.gz" \
	-G StandardAnnotation \
	--use-new-qual-calculator \
	--tmp-dir=${WORKDIR}
}
