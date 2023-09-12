#!/bin/bash
## This script runs the joint pipeline on a single interval
# The three arguments are the COHORT, the number of intervals, and the interval to be run
# The COHORT argument is the name of the location in /storage1/fs1/${STORAGE_USER}/Active/$USER/c1in
# The NUM_INTERVALS argument is the number of intervals in the run
# The INTERVAL is the interval to run because it failed elsewhere

COHORT=$1
NUM_INTERVALS=$2
INTERVAL=$3

# These variables can be changed to run for other users
export COMPUTE_USER=fernandezv
export STORAGE_USER=cruchagac
export SCRATCH_USER=cruchagac
REF_DIR=/scratch1/fs1/cruchagac/WXSref

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_FILE=$(bash ${SCRIPT_DIR}/..emacs/makeCohortEnv5.bash $COHORT $NUM_INTERVALS)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

## 1a. Joint Call on Intervals - genomicsDB and GenotypeGVCFs
# This job runs a single interval that failed
# This job adds the intervals from all input gvcfs to a genomicsdb.
# After the genomicsDBImport is completed, the job calls variants with GenotypeGVCFs
# This job results in a single vcf for the interval given
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-call-rescue-$INTERVAL \
    -Ne \
    -sp 70 \
    -n 1 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.${INTERVAL}.joint_s1.%J.out \
    -R 'select[mem>80GB && tmp>10GB] rusage[mem=80GB,tmp=10GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointcaller:2.0)' \
    bash /scripts/jointcallsplitinterval.bash $INTERVAL
