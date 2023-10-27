#!/bin/bash
## This script runs the joint pipeline on a single chromosome
# The two arguments are the COHORT and the CHR to be run
# The COHORT argument is the name of the location in /storage1/fs1/${STORAGE_USER}/Active/$USER/c1out
# The CHR argument is the chromosome to be run (ex. chr1, chrX)
COHORT=$1
CHR=$2

# These variables can be changed to run for other users
export COMPUTE_USER=fernandezv
export STORAGE_USER=cruchagac
export SCRATCH_USER=cruchagac
REF_DIR=/scratch1/fs1/cruchagac/WXSref

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_FILE=$(bash ${SCRIPT_DIR}/../makeCohortEnv.bash $COHORT $CHR)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}-${CHR}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

## 3. Joint QC
# This step performs VQSR filtering and a host of other filters on the CHR joint vcf
# The outcome from this step is a file containing counts of the variants in each step and a final filtered joing vcf for the CHR
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-qc \
    -N \
    -n 4 \
    -sp 90 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.${CHR}.joint_s3.%J.out \
    -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointasqc:2.0)' \
    bash /scripts/VQCPipelinekeep.bash
