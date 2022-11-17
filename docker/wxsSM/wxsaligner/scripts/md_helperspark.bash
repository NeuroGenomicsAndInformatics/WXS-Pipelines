#!/bin/bash
THREADS=$(( LSB_MAX_NUM_PROCESSORS * 2 ))
MD_INPUTS=()
for BM in $(find $INDIR -name "*.f*q*.bam"); do
ln -s $BM /tmp/${BM##*/}
MD_INPUTS+="-I /tmp/${BM##*/} "
done
${GATK} \
  --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  MarkDuplicatesSpark \
    ${MD_INPUTS[@]}\
    -O ${OUTDIR}/marked.bam \
    --conf "spark.executor.cores=$(( THREADS - 2 ))" \
    --conf "spark.local.dir=$TMP_DIR" \
&& rm -R ${INDIR} \
&& ${GATK} \
  --java-options "-Xmx170g -XX:ParallelGCThreads=2 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SetNmMdAndUqTags \
    -I ${OUTDIR}/marked.bam \
    -O ${OUTDIR}/${CRAM} \
    -R ${REF_FASTA} \
    --SET_ONLY_UQ true
rm -R ${INDIR}
rm ${OUTDIR}/marked.bam
samtools index -@ $LSB_MAX_NUM_PROCESSORS $OUTDIR/$CRAM
