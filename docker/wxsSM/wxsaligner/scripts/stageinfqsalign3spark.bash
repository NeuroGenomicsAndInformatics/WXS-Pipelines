#!/bin/bash
rsync -rL $STAGE_INDIR/ $INDIR
INFQ_FILE=${INDIR}/infqfile.txt
echo -n "" > $INFQ_FILE
for FQ in $(find $INDIR -name "*_1.fastq.gz"); do
SM=$(echo $FULLSMID | cut -d^ -f1)
BARCODE=$(echo $FULLSMID | cut -d^ -f2)
PROJECT=$(echo $FULLSMID | cut -d^ -f3)
FLOWCELL=$(echo ${FQ##*/} | cut -d_ -f1 | cut -d. -f1)
LANE=$(echo ${FQ##*/} | cut -d_ -f1 | cut -d. -f2)
echo "@RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" > ${OUTDIR}/${FULLSMID}.${FLOWCELL}_${LANE}.rgfile
echo "${FQ} ${FQ/_1.fastq/_2.fastq} @RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >> ${INFQ_FILE}
done
JOBS_IN_ARRAY=$(wc -l ${INDIR}/infqfile.txt | cut -d ' ' -f1)
LSF_DOCKER_VOLUMES="/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -q general \
  -g /matthew.j/compute-fernandezv \
  -G compute-fernandezv \
  -J ngi-$USER-$FULLSMID-fqsalign[1-${JOBS_IN_ARRAY}] \
  -n 4 \
  -sp 55 \
  -o ${LOGDIR}/align_${FULLSMID}.%J.%I.out \
  -R 'select[mem>80GB] rusage[mem=80GB] span[hosts=1]' \
  -a 'docker(mjohnsonngi/wxsalignhelper:2.0)' \
  bash /scripts/bwa_helperfqs4spark.bash
bash bwa_helperfqs4sparkbase.bash
bwait -w "done(\"ngi-$USER-$FULLSMID-fqsalign\")" \
&& bash md_helperspark.bash
