#!/bin/bash
export THREADS=8
export MEM=96
JOB_GROUP="/${USER}/compute-cruchagac"
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/matthewj/baseEnvs/pipelinebase.env /scratch1/fs1/cruchagac/matthewj/baseEnvs/references.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-test \
-n ${THREADS} \
-o /scratch1/fs1/cruchagac/matthewj/c1out/DEBUG_TEST/test_s1.%J.out \
-e /scratch1/fs1/cruchagac/matthewj/c1out/DEBUG_TEST/test_s1.%J.err \
-R 'select[mem>102000 && tmp>10] rusage[mem=100000, tmp=10] span[hosts=1]' \
-M 120000000 \
-G compute-cruchagac \
-q general-interactive \
-a 'docker(mjohnsonngi/pipelinea:stable)' /bin/bash "$1"
