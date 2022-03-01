#!/bin/bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac/${USER}/c1in:/input \
/scratch1/fs1/cruchagac/WXSref:/ref \
/scratch1/fs1/cruchagac/${USER}/c1out:/output \
/storage1/fs1/cruchagac/Active/${USER}/c1out:/final_output"
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
for FULLSMID in $(cat $1); do
  bash ./makeSampleEnv.bash $FULLSMID
  LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/${USER}/c1in/envs/${FULLSMID}.env" \
  bsub -g ${JOB_GROUP} \
  -J ngi-${USER}-stage2-$FULLSMID \
  -n 8 \
  -N \
  -o /scratch1/fs1/cruchagac/${USER}/c1out/logs/${FULLSMID}/${FULLSMID}_s2.%J.out \
  -R 'select[mem>105000] rusage[mem=105000] span[hosts=1]' \
  -M 110000 \
  -G compute-cruchagac \
  -q general \
  -a 'docker(mjohnsonngi/wxspipeline:dev)' /scripts/pipelineStage2.bash
  done
