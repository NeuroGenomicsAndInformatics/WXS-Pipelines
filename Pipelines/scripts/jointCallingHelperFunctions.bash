function stageDataForCOHORT ()
{
	rsync -rL $STAGE_INDIR/ $INDIR
	for FILE in $(ls $INDIR); do
		[[ ! -h ${INDIR}/${FILE//^/.} ]] && ln -s ${INDIR}/${FILE} ${INDIR}/${FILE//^/.}
	done
}
function transferOutputFilesToStorage ()
{
	[ ! -d ${FINAL_OUTDIR} ] && mkdir ${FINAL_OUTDIR}
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
	SAMPLE_MAP="${WORKDIR}/SampleMap.txt"
	echo -n "" > ${SAMPLE_MAP}
	for SAMPLE in $(find ${INDIR} -name "*.vcf.gz" ${pwd} | grep "\^") ; do
		echo -e "$(echo ${SAMPLE##*/} | cut -d'.' -f1)\tfile://${SAMPLE//^/.}" >> ${SAMPLE_MAP}
	done
}
function buildGenomicDB ()
{
	DATABASE="${WORKDIR}/db${INTERVAL}"
[ ! -d ${OUTDIR}/dbwork ] && mkdir ${OUTDIR}/dbwork
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S3MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
	GenomicsDBImport \
	-L ${INTERVAL} \
	--sample-name-map ${SAMPLE_MAP} \
	--genomicsdb-workspace-path ${WORKDIR}/db${INTERVAL} \
	--batch-size 50 \
	--tmp-dir $OUTDIR/dbwork \
	--reader-threads ${S3THREADS}
#	--batch-size 100 \
#	--merge-input-intervals true \
#	--genomicsdb-shared-posixfs-optimizations true \
}
function jointCallCohort ()
{
	JOINT_GVCF="${OUTDIR}/${COHORT}.${INTERVAL}.joint.g.vcf.gz"
${GATK} --java-options "-Xms${MEM_SPLIT}g -Xmx${S3MEM}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
GenotypeGVCFs \
	-R ${REF_FASTA} \
	-V gendb://${DATABASE} \
	-O "${WORKDIR}/${COHORT}.${INTERVAL}.joint.g.vcf.gz" \
	-G StandardAnnotation \
	--tmp-dir ${WORKDIR}
#	--genomicsdb-shared-posixfs-optimizations true
rsync ${WORKDIR}/${COHORT}.${INTERVAL}.joint.g.vcf.gz* ${OUTDIR}/
}
