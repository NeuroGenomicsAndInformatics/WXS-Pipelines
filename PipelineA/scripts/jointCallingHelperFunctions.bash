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
function cleanup ()
{
	#rm -R $STAGE_INDIR
	rm -R $INDIR
	#rm -R $OUTDIR
}
function makeSampleMap ()
{
	SAMPLE_MAP="${WORKDIR}/SampleMap.txt"
	touch ${SAMPLE_MAP}
	for SAMPLE in $(ls -d ${INDIR}) ; do
	echo -en "${echo $SAMPLE | cut -d'^' -f1}\t" >> ${SAMPLE_MAP}
	echo "$(find ${STAGE_DIR}/${SAMPLE} -name "*.vcf.gz" -print)" >> ${SAMPLE_MAP}
	done
}
function buildGenomicDB ()
{
	DATABASE="${WORKDiR}/db"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	GenomicDBImport \
	--sample-name-map ${SAMPLE_MAP} \
	--genomicsdb-workspace-path ${WORKDIR}/db \
	--batch-size 50 \
	--tmp-dir $OUTDIR/dbwork \
	--reader-threads ${S3THREADS}
}
function jointCallCohort ()
{
	JOINT_GVCF="${OUTDIR}/${COHORT}.joint.g.vcf.gz"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
GenotypeGVCFs \
	-R ${REF_FASTA} \
	-V gendb://${DATABASE} \
	-O "${OUTDIR}/${COHORT}.joint.g.vcf.gz" \
	-G StandardAnnotation \
	--use-new-qual-calculator \
	--tmp-dir=${WORKDIR}
}
