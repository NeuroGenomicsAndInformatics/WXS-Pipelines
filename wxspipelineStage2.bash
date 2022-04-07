#!/bin/bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
  /scratch1/fs1/cruchagac/${USER}/c1in:/input \
  /scratch1/fs1/cruchagac/WXSref:/ref \
  /scratch1/fs1/cruchagac/${USER}/c1out:/output \
  /storage1/fs1/cruchagac/Active/${USER}/c1out:/final_output"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
for FULLSMID in $(cat $1); do
  bash ${SCRIPT_DIR}/makeSampleEnv.bash $FULLSMID
  LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/${USER}/c1in/envs/${FULLSMID}.env" \
  bsub -g ${JOB_GROUP} \
  -J ngi-${USER}-stage2-${FULLSMID} \
  -N \
  -o /scratch1/fs1/cruchagac/${USER}/c1out/logs/${FULLSMID}/${FULLSMID}_s2.%J.out \
  -R 'select[mem>105GB && ncpus>8] rusage[mem=105GB] span[hosts=1]' \
  -M 110GB \
  -G compute-cruchagac \
  -q general \
  -a 'docker(mjohnsonngi/wxspipeline:1.1)' /scripts/pipelineStage2.bash
  done
