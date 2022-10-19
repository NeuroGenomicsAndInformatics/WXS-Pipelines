#!/bin/bash
TMPDIR=/scratch1/fs1/cruchagac/$USER/c1out/${FULLSMID//^/}
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
  -n 8 \
  -o ${LOGDIR}/align_${FULLSMID}.%J.%I.out \
  -R 'select[mem>160GB] rusage[mem=160GB] span[hosts=1]' \
  -a 'docker(mjohnsonngi/wxsalignhelper:2.0)' \
  bash /scripts/bwa_helperfqs4.bash
LSF_DOCKER_VOLUMES="/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
/scratch1/fs1/cruchagac/WXSref:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g /matthew.j/compute-fernandezv \
  -J ngi-$USER-$FULLSMID-md \
  -w "done(\"ngi-$USER-$FULLSMID-fqsalign\")" \
  -n 4 \
  -K \
  -sp 70 \
  -o ${LOGDIR}/md.%J.out \
  -R 'select[mem>240GB] rusage[mem=240GB] span[hosts=1]' \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxsalignhelper:2.0)' \
  bash /scripts/md_helper.bash \
&& echo "$FULLSMID aligned, sorted, and marked."
