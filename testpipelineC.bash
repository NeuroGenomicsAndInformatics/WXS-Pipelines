#!/bin/bash
export THREADS=8
export MEM=32
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
cat $1 | while read LINE; do
CRAM="$(echo $LINE | cut -d ' ' -f2)"
bash ${SCRIPT_DIR}/makeCramEnv.bash $LINE
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac/${USER}/c1in:/input \
/scratch1/fs1/cruchagac/WXSref:/ref \
/scratch1/fs1/cruchagac/${USER}/c1out:/output \
/storage1/fs1/cruchagac/Active/${USER}/c1out:/final_output" \
LSF_DOCKER_ENV_FILE="/scratch1/fs1/cruchagac/${USER}/c1in/envs/${CRAM}.env" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage0-$CRAM \
-n ${THREADS} \
-cwd ${SCRIPT_DIR} \
-o /scratch1/fs1/cruchagac/${USER}/c1out/${CRAM%.cram}/${CRAM%.cram}.%J.out \
-R 'select[mem>32000] rusage[mem=32000]' \
-M 46000000 \
-G compute-cruchagac \
-q general \
-a 'docker(mjohnsonngi/wxspipeline:dev)' bash /scripts/pipelineCStage0.bash $CRAM
done
