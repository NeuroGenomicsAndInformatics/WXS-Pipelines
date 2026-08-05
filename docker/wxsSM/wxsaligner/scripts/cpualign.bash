#!/bin/bash
JOBS_IN_ARRAY=$(wc -l ${INDIR}/infqfile.txt | cut -d ' ' -f1)
for (( i=0 ; $i < $JOBS_IN_ARRAY ; i++ )); do
  touch ${INDIR}/${i}.lock
done
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -q general \
  -g /matthew.j/compute-${COMPUTE_USER} \
  -G compute-${COMPUTE_USER} \
  -J fqsalign${LSB_JOBID}[1-$(( JOBS_IN_ARRAY - 1 ))] \
  -n 4 \
  -sp 55 \
  -o ${LOGDIR}/align_${FULLSMID}.%J.%I.out \
  -R 'select[mem>120GB] rusage[mem=120GB] span[hosts=1]' \
  -a 'docker(mjohnsonngi/wxsalignhelper:2.1)' \
  bash /scripts/bwa_helperfqs4.bash
bash /scripts/bwa_helperfqs4base.bash
bsub -q general \
  -g /matthew.j/compute-${COMPUTE_USER} \
  -G compute-${COMPUTE_USER} \
  -w "done(fqsalign${LSB_JOBID}*)" -ti \
  -n 1 \
  -Ne \
  -K \
  -sp 55 \
  -a 'docker(mjohnsonngi/wxsalignhelper:2.1)' \
  'echo "hello" > ${INDIR}/holder_exit.txt'
[[ -z $(find ${INDIR} -maxdepth 1 -name "*.lock") ]] || [[ -f /${INDIR}/holder_exit.txt ]] || exit 67
bash /scripts/md_helper.bash