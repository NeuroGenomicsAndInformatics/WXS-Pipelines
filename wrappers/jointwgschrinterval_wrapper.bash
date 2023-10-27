#!/bin/bash
## This script runs the joint calling portion of the pipeline on a single interval
# The script should be used when one of the jobs in an array exits
# The three arguments are the COHORT, the CHR, and the INTERVAL to be run
# The COHORT argument is the name of the location in /storage1/fs1/${STORAGE_USER}/Active/$USER/c1in
# The CHR argument is the chromosome to be run (ex. chr1, chrX)
# The INTERVAL argument is the job number in the array that failed
COHORT=$1
CHR=$2
INTERVAL=$3

# These variables can be changed to run for other users
export COMPUTE_USER=fernandezv
export STORAGE_USER=cruchagac
export SCRATCH_USER=cruchagac
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_FILE=$(bash ${SCRIPT_DIR}/../makeCohortEnv.bash $COHORT $CHR)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}-${CHR}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint-rescue"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

# This component is to find the number of jobs to use in an array for calling
SHARDS=50
echo $SHARDS

## 1a. Joint Call on Intervals - genomicsDBImport & GenotypeGVCFs
# This first job is an array of jobs based on the number of intervals in the interval list for the CHR given
# These jobs each add the intervals from all input gvcfs to a genomicsdb
# This genomicsdb is then used to joint call the variants in the interval
# The outcome of this step is a joint vcf the covers a single interval
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-call-rescue-$INTERVAL \
    -N \
    -sp 90 \
    -n 2 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.${CHR}.joint_s1.%J.%I.out \
    -R 'select[mem>120GB] rusage[mem=120GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointcaller:2.0)' \
    bash /scripts/jointcallchrsplitinterval_rescue.bash $INTERVAL
