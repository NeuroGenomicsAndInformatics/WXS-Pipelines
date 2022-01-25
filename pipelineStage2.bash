#!/bin/bash
export THREADS=8
export MEM=48
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
bash ./perSampleEnvs.bash $1
for FULLSMID in $(cat $1); do
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/baseEnvs/pipelineBase.env /scratch1/fs1/cruchagac/matthewj/baseEnvs/references.env /scratch1/fs1/cruchagac/matthewj/c1in/envs/${FULLSMID}.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage2-$FULLSMID \
-n ${THREADS} \
-R 'select[mem>48000] rusage[mem=48GB]' \
-M ${MEM}GB \
-W 720 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/pipelinea:latest)' /scripts/pipelineABStage2.bash
done
