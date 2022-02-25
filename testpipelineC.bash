#!/bin/bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active/${USER}/c1in:/staged_input \
/scratch1/fs1/cruchagac/${USER}/c1in:/input \
/scratch1/fs1/cruchagac/WXSref:/ref \
/scratch1/fs1/cruchagac/${USER}/c1out:/output \
/storage1/fs1/cruchagac/Active/${USER}/c1out:/final_output"
export THREADS=8
export MEM=32
export CRAM="$1"
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
LSF_DOCKER_ENV_FILE="./baseEnvs/pipelinebase.env ./baseEnvs/references.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage0-$CRAM \
-n ${THREADS} \
-o /scratch1/fs1/cruchagac/${USER}/c1out/${CRAM%.cram}/${CRAM%.cram}.%J.%I.out \
-e /scratch1/fs1/cruchagac/${USER}/c1out/${CRAM%.cram}/${CRAM%.cram}.%J.%I.err \
-R 'select[mem>32000] rusage[mem=32000] span[hosts=1]' \
-M 46000000 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/pipelinea:latest)' /scripts/pipelineCStage0.bash
