#!/bin/bash
export LSF_DOCKER_VOLUMES="/scratch1/fs1/cruchagac/matthewj/c1in:/input \
/scratch1/fs1/cruchagac/matthewj/ref:/ref \
/scratch1/fs1/cruchagac/matthewj/c1out:/output \
/storage1/fs1/cruchagac/Active/matthewj/c1out:/final_output"
export THREADS=16
export MEM=96
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
bash ./perSampleEnvs.bash $1
for FULLSMID in $(cat $1); do
JOBS_IN_ARRAY=$(ls /scratch1/fs1/cruchagac/matthewj/c1in/${FULLSMID}/*.rgfile | wc -w)
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/c1in/envs/${FULLSMID}.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage1-$FULLSMID[1-$JOBS_IN_ARRAY] \
-n ${THREADS} \
-o /scratch1/fs1/cruchagac/matthewj/c1out/${FULLSMID}/${FULLSMID}_s1.%J.%I.out \
-e /scratch1/fs1/cruchagac/matthewj/c1out/${FULLSMID}/${FULLSMID}_s1.%J.%I.err \
-R 'select[mem>102000 && tmp>10] rusage[mem=100000, tmp=10]' \
-M 120000000 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/wxspipeline:latest)' /scripts/pipelineStage1.bash
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/c1in/envs/${FULLSMID}.env" \
bsub -g ${JOB_GROUP} \
-w "done(\"ngi-${USER}-stage1-$FULLSMID\")" \
-J ngi-${USER}-stage2-$FULLSMID \
-n ${THREADS} \
-N \
-o /scratch1/fs1/cruchagac/matthewj/c1out/${FULLSMID}/${FULLSMID}_s2.%J.out \
-e /scratch1/fs1/cruchagac/matthewj/c1out/${FULLSMID}/${FULLSMID}_s2.%J.err \
-R 'select[mem>102000] rusage[mem=100000] span[hosts=1]' \
-M 120000000 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/wxspipeline:latest)' /scripts/pipelineStage2.bash
done
