#!/bin/bash
export THREADS=8
export MEM=48
JOB_GROUP="/${USER}/compute-cruchagac"
FULLSM="$(head -n 1 $1)"
bash ./perSampleEnvs.bash $1
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/c1in/envs/${FULLSM}.env" \
bsub -g ${JOB_GROUP} -Is \
-J ngi-${USER}-test \
-n ${THREADS} \
-R 'select[mem>48000] rusage[mem=48GB]' \
-M ${MEM}GB \
-G compute-cruchagac \
-q general-interactive \
-a 'docker(mjohnsonngi/pipelinea:stable)' /bin/bash
