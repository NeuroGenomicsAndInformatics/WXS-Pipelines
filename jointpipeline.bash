#!/bin/bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active/$USER/c1in:/staged_input \
/scratch1/fs1/cruchagac/$USER/c1in:/input \
/scratch1/fs1/cruchagac/$USER/ref:/ref \
/scratch1/fs1/cruchagac/$USER/c1out:/output \
/storage1/fs1/cruchagac/Active/$USER/c1out:/final_output"
JOB_GROUP="/${USER}/compute-cruchagac"
COHORT="$1"
bgadd -L 10 ${JOB_GROUP}
bash ./makeCohortEnv.bash ${COHORT} $2
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/c1in/envs/${FULLSMID}.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage3-$COHORT \
-n 16 \
-o /scratch1/fs1/cruchagac/$USER/c1out/${COHORT}/${COHORT}_s1.%J.%I.out \
-R 'select[mem>150000] rusage[mem=150000/job] span[hosts=1]' \
-M 160000 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/wxspipeline:joint)' /scripts/jointCalling.bash
