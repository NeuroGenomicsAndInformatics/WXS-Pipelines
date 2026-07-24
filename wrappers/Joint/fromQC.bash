#!/bin/bash
## This script runs the joint pipeline on a the entire genome by splitting it into intervals
# The two arguments are the COHORT and the NUM_INTERVALS to split the genome into
# The COHORT argument is the name of the location in /storage1/fs1/${STORAGE_USER}/Active/$USER/c1in
# The NUM_INTERVALS argument is the number of intervals to split the genome into
COHORT=$1
NUM_INTERVALS=$2

# These variables can be changed to run for other users
export COMPUTE_USER=fernandezv
export STORAGE_USER=cruchagac
export SCRATCH_USER=cruchagac
REF_DIR=/scratch1/fs1/cruchagac/WXSref

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_FILE=$(bash ${SCRIPT_DIR}/../../makeCohortEnvInt.bash $COHORT $NUM_INTERVALS)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

## 5. Gather QCed Vcfs
# This step takes all of the joint vcfs from the previous step and combines them into a chromosome joint vcf
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-gather-qced \
    -N \
    -n 4 \
    -sp 80 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s5.%J.out \
    -R 'select[mem>220GB] rusage[mem=220GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointgatherer:2.1)' \
    bash /scripts/gathervcfs_called.bash
