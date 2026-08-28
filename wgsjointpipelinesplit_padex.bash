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
ENV_FILE=$(bash ${SCRIPT_DIR}/makeCohortEnvInt.bash $COHORT $NUM_INTERVALS)

# Pipeline variable setup for running the jobs
JOBNAME="ngi-${USER}-${COHORT}"
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}/joint"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 50 ${JOB_GROUP}
[ ! -d /scratch1/fs1/${COMPUTE_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

# The first calling stage is looped to run all of the splits
# This script doesn't use job arrays because of the max of 1000 jobs in an array
for ((i=0 ; i < $NUM_INTERVALS ; i++)); do

## 1. Joint Call on Intervals - genomicsDB and GenotypeGVCFs
# This first job is split into the number of intervals with each interval having its own job
# These jobs each add the intervals from all input gvcfs to a genomicsdb.
# After the genomicsDBImport is completed, the job calls variants with GenotypeGVCFs
# This job results in a single vcf for the interval given
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-call-$i \
    -Ne \
    -sp 70 \
    -n 1 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s1.%J.${i}.out \
    -R 'select[mem>80GB && tmp>6GB] rusage[mem=80GB,tmp=6GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointcaller:2.1)' \
    bash /scripts/jointcallsplitinterval_padex.bash $i
done

## 2. Gather Joint Called Vcfs
# This step takes all of the joint vcfs from the previous step and combines them into a chromosome joint vcf
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -w "done(${JOBNAME}-call-*)" -ti \
    -J ${JOBNAME}-gather-called \
    -N \
    -n 4 \
    -sp 80 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s2.%J.out \
    -R 'select[mem>220GB] rusage[mem=220GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointgatherer:2.1)' \
    bash /scripts/gathervcfs_called.bash

## 3. VQSR
# This step performs VQSR calculations on the joint called dataset
# The outcome from this step are 2 recal table files. One for SNP and one for INDEL.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -w "done(${JOBNAME}-gather-called)" -ti \
    -J ${JOBNAME}-vqsr-snp \
    -N \
    -n 4 \
    -sp 90 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s3_snp.%J.out \
    -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointasqc:2.1)' \
    bash /scripts/VQSR_SNP.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -w "done(${JOBNAME}-gather-called)" -ti \
    -J ${JOBNAME}-vqsr-indel \
    -N \
    -n 4 \
    -sp 90 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s3_indel.%J.out \
    -R 'select[mem>100GB] rusage[mem=100GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointasqc:2.1)' \
    bash /scripts/VQSR_INDEL.bash

## 4. QC on Intervals
# These jobs QC the joint called intervals.
# They use the split vcfs from the calling stage.
# First it applies the VQSR done on the gathered vcf.
# It then does several hard filtering steps.

for ((i=0 ; i < $NUM_INTERVALS ; i++)); do

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-qc-$i \
    -w "done(\"${JOBNAME}-vqsr*\")" -ti \
    -Ne \
    -sp 70 \
    -n 1 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s4.%J.${i}.out \
    -R 'select[mem>80GB && tmp>6GB] rusage[mem=80GB,tmp=6GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointasqc:2.1)' \
    bash /scripts/VQCPipeline.bash $i

done

## 5. Gather QCed Vcfs
# This step takes all of the joint vcfs from the previous step and combines them into a chromosome joint vcf
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
    -w "done(${JOBNAME}-qc-*)" -ti \
    -J ${JOBNAME}-gather-called \
    -N \
    -n 4 \
    -sp 80 \
    -o /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}.joint_s5.%J.out \
    -R 'select[mem>220GB] rusage[mem=220GB] span[hosts=1]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointgatherer:2.1)' \
    bash /scripts/gathervcfs_qced.bash
