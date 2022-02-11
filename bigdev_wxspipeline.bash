#!/bin/bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active/matthewj/c1in/staged_input:/staged_input \
/scratch1/fs1/cruchagac/matthewj/c1in:/input \
/scratch1/fs1/cruchagac/matthewj/ref:/ref \
/scratch1/fs1/cruchagac/matthewj/c1out:/output \
/storage1/fs1/cruchagac/Active/matthewj/c1out:/final_output"
export THREADS=16
export MEM=96
JOB_GROUP="/${USER}-dev/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
bash ./perSampleEnvs.bash $1
for FULLSMID in $(cat $1); do
JOBS_IN_ARRAY=$(ls /scratch1/fs1/cruchagac/matthewj/c1in/${FULLSMID}/*.rgfile | wc -w)
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/c1in/envs/${FULLSMID}.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage1-$FULLSMID[1-$JOBS_IN_ARRAY] \
-n 4 \
-o /scratch1/fs1/cruchagac/matthewj/c1out/${FULLSMID}/${FULLSMID}_s1.%J.%I.out \
-e /scratch1/fs1/cruchagac/matthewj/c1out/${FULLSMID}/${FULLSMID}_s1.%J.%I.err \
-R 'select[mem>8000] rusage[mem=8000]' \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/wxspipeline:dev)' /scripts/bigtest_pipelineStage1.bash
done
