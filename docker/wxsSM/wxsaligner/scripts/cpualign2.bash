#!/bin/bash
JOBS_IN_ARRAY=$(wc -l ${INDIR}/infqfile.txt | cut -d ' ' -f1)
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -q general \
  -g /matthew.j/compute-fernandezv \
  -G compute-fernandezv \
  -J ngi-$USER-$FULLSMID-fqsalign[1-$(( JOBS_IN_ARRAY ))] \
  -n 4 \
  -sp 55 \
  -o ${LOGDIR}/align_${FULLSMID}.%J.%I.out \
  -R 'select[mem>90GB] rusage[mem=90GB] span[hosts=1]' \
  -a 'docker(mjohnsonngi/wxsalignhelper:2.0)' \
  bash /scripts/bwa_helperfqs4.bash
bwait -w "done(\"ngi-$USER-$FULLSMID-fqsalign\")" || exit 67
bash /scripts/md_helper.bash
