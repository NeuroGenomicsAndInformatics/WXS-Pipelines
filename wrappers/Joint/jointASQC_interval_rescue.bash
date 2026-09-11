#!/bin/bash
## This script runs the QC pipeline on a single vcf
# The three arguments are the COHORT, the NUM_INTERVALS, and INTERVAL.
# The COHORT argument is the name of the location in /storage1/fs1/${STORAGE_USER}/Active/$USER/c1out
# The NUM_INTERVALS argument is the number of intervals in the joint call.
# The INTERVAL is the interval to be run.
COHORT=$1
NUM_INTERVALS=$2
INTERVAL=$3

# These variables can be changed to run for other users
export COMPUTE_USER=fernandezv
export STORAGE_USER=cruchagac
export SCRATCH_USER=cruchagac
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_FILE=$(bash ${SCRIPT_DIR}/../../makeCohortEnvInt.bash $COHORT $NUM_INTERVALS)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}-${CHR}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

## 4. Joint QC
# This step performs VQSR filtering and a host of other filters on the CHR joint vcf
# The outcome from this step is a file containing counts of the variants in each step and a final filtered joing vcf for the CHR
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-qc \
    -Ne \
    -n 4 \
    -sp 90 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s4.${INTERVAL}.%J.out \
    -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointasqc:2.1)' \
    bash /scripts/VQCPipeline.bash $INTERVAL
