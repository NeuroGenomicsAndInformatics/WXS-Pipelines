#!/bin/bash
# Check if running in the background. If not, call itself in the background and exit the foreground version
export THREADS=8
export MEM=32
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
JOBS_IN_ARRAY=$(cat $1 | wc -l)
for CRAM in $(cat $1); do
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
bwait -w "done(\"ngi-${USER}-stage0-$CRAM\")" \
&& bash ./pipelineA.bash ~/work$CRAM.txt &
done
disown
