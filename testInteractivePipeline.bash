#!/bin/bash
export THREADS=16
export MEM=48
JOB_GROUP="/${USER}/compute-cruchagac"
FULLSM="$(head $1)"
bash ./perSampleEnvs.bash $1
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/baseEnvs/pipelineBase.env /scratch1/fs1/cruchagac/matthewj/baseEnvs/references.env /scratch1/fs1/cruchagac/matthewj/c1in/envs/${FULLSM}.env" \
bsub -g ${JOB_GROUP} -Is \
-J ngi-${USER}-test \
-n ${THREADS} \
-R 'select[mem>48000] rusage[mem=48GB]' \
-M ${MEM}GB \
-G compute-cruchagac \
-q general-interactive \
-a 'docker(mjohnsonngi/pipelinea:latest)' /bin/bash
