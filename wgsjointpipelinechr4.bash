#!/bin/bash
## This script runs the joint pipeline on a single chromosome
# The two arguments are the COHORT and the CHR to be run
# The COHORT argument is the name of the location in /storage1/fs1/${STORAGE_USER}/Active/$USER/c1in
# The CHR argument is the chromosome to be run (ex. chr1, chrX)
COHORT=$1
CHR=$2

# These variables can be changed to run for other users
export COMPUTE_USER=fernandezv
export STORAGE_USER=cruchagac
export SCRATCH_USER=cruchagac
REF_DIR=/scratch1/fs1/fernandezv/WXSref

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENV_FILE=$(bash ${SCRIPT_DIR}/makeCohortEnv.bash $COHORT $CHR)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}-${CHR}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

# This component is to find the number of jobs to use in an array for calling
# The number echoed represents the number of intervals in the interval list for that CHR
SHARDS=50
echo $SHARDS

## 1a. Joint Call on Intervals - genomicsDB
# This first job is an array of jobs based on the number of intervals in the interval list for the CHR given
# These jobs each add the intervals from all input gvcfs to a genomicsdb.
# The outcome from this step should be a genomicsDB for each interval
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-call[1-${SHARDS}] \
    -Ne \
    -sp 70 \
    -n 2 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.${CHR}.joint_s1.%J.%I.out \
    -R 'select[mem>240GB] rusage[mem=240GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointcaller:2.0)' \
    bash /scripts/jointcallinterval4.bash

## 2. Sort Vcfs
# This step takes all of the joint vcfs from the previous step and combines them into a chromosome joint vcf
# The sorting may be mostly unnecessary, but it does eliminate a problem if the interval vcfs are out of order
# The outcome of this step is a single joint vcf that covers the CHR given
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -w "done(${JOBNAME}-call)" \
    -J ${JOBNAME}-sort \
    -N \
    -n 4 \
    -sp 80 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.${CHR}.joint_s2.%J.out \
    -R 'select[mem>220GB] rusage[mem=220GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointsorter:2.0)' \
    bash /scripts/mergevcfs.bash

## 3. Joint QC
# This step performs VQSR filtering and a host of other filters on the CHR joint vcf
# The outcome from this step is a file containing counts of the variants in each step and a final filtered joing vcf for the CHR
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -w "done(${JOBNAME}-sort)" \
    -J ${JOBNAME}-qc \
    -N \
    -n 4 \
    -sp 90 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.${CHR}.joint_s3.%J.out \
    -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointasqc:2.0)' \
    bash /scripts/VQCPipeline.bash
