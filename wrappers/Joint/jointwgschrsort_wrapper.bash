#!/bin/bash
## This script runs the sorting portion of the pipeline on a single chromosome
# The script should be used when one of the jobs in an array exits
# The two arguments are the COHORT and the CHR
# The COHORT argument is the name of the location in /storage1/fs1/${STORAGE_USER}/Active/$USER/c1in
# The CHR argument is the chromosome to be run (ex. chr1, chrX)
COHORT=$1
CHR=$2

# These variables can be changed to run for other users
export COMPUTE_USER=fernandezv
export STORAGE_USER=cruchagac
export SCRATCH_USER=cruchagac
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_FILE=$(bash ${SCRIPT_DIR}/../../makeCohortEnv.bash $COHORT $CHR)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}-${CHR}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint-rescue"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

## 2. Sort Vcfs
# This step takes all of the joint vcfs from the previous step and combines them into a chromosome joint vcf
# The sorting may be mostly unnecessary, but it does eliminate a problem if the interval vcfs are out of order
# The outcome of this step is a single joint vcf that covers the CHR given
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-sort-rescue \
    -N \
    -n 4 \
    -sp 80 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.${CHR}.joint_s2.%J.out \
    -R 'select[mem>220GB] rusage[mem=220GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointsorter:2.0)' \
    bash /scripts/sortvcfs.bash
