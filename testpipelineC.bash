#!/bin/bash
# Check if running in the background. If not, call itself in the background and exit the foreground version
export THREADS=1
export MEM=32
export CRAM="$1"
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
LSF_DOCKER_ENV_FILE="./baseEnvs/pipelinebase.env ./baseEnvs/references.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage0-$CRAM \
-n ${THREADS} \
-o /scratch1/fs1/cruchagac/matthewj/c1out/${CRAM%.cram}/${CRAM%.cram}.%J.%I.out \
-e /scratch1/fs1/cruchagac/matthewj/c1out/${CRAM%.cram}/${CRAM%.cram}.%J.%I.err \
-R 'select[mem>32000] rusage[mem=32000] span[hosts=1]' \
-M 46000000 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/pipelinea:latest)' /scripts/pipelineCStage0.bash
done
