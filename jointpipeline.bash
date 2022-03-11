#!/bin/bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac/$USER/c1in:/input \
/scratch1/fs1/cruchagac/WXSref:/ref \
/scratch1/fs1/cruchagac/$USER/c1out:/output \
/storage1/fs1/cruchagac/Active/$USER/c1out:/final_output"
JOB_GROUP="/${USER}/compute-cruchagac"
COHORT="$1"
INTERVAL="$2"
bgadd -L 10 ${JOB_GROUP}
if [[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/logs ]]; then mkdir /scratch1/fs1/cruchagac/${USER}/c1out/logs; fi
bash ./makeCohortEnv.bash ${COHORT} ${INTERVAL}
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/${USER}/c1in/envs/${COHORT}_${INTERVAL}.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage3-${COHORT}_${INTERVAL} \
-N \
-n 4 \
-sp 55 \
-o /scratch1/fs1/cruchagac/${USER}/c1out/logs/${COHORT}_${INTERVAL}/${COHORT}_${INTERVAL}_s1.%J.%I.out \
-R 'select[mem>80000 && tmp>100] rusage[mem=80000/job, tmp=100] span[hosts=1]' \
-M 82000 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/wxspipeline:joint)' /scripts/jointCalling.bash
