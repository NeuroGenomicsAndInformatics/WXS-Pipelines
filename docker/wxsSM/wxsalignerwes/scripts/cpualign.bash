#!/bin/bash
JOBS_IN_ARRAY=$(wc -l ${INDIR}/infqfile.txt | cut -d ' ' -f1)
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -q general \
  -g /matthew.j/compute-${COMPUTE_USER} \
  -G compute-${COMPUTE_USER} \
  -J ngi-$USER-$FULLSMID-fqsalign[1-$(( JOBS_IN_ARRAY - 1 ))] \
  -n 4 \
  -sp 55 \
  -o ${LOGDIR}/align_${FULLSMID}.%J.%I.out \
  -R 'select[mem>80GB] rusage[mem=80GB] span[hosts=1]' \
  -a 'docker(mjohnsonngi/wxsalignhelper:2.0)' \
  bash /scripts/bwa_helperfqs4wes.bash
bash /scripts/bwa_helperfqs4base.bash
[[ $JOBS_IN_ARRAY -eq 1 ]] || bwait -w "done(\"ngi-$USER-$FULLSMID-fqsalign\")" || exit 67
bash /scripts/md_helper.bash
