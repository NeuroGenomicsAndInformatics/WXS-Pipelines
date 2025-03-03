#!/bin/bash
## Set up variables for specific users. These should be all that's needed to change user
export COMPUTE_USER=fernandezv
export SCRATCH_USER=cruchagac
export STORAGE_USER=cruchagac
export REF_DIR="/scratch1/fs1/cruchagac/WXSref"

## 0. Set up for job submission 
# 0.1 Make expected directories in case they are missing
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

# 0.2 Priorities are set to handle bounded-buffer issues
PRIORITY_INTLIST=55
PRIORITY_MISS=60
PRIORITY_ANNAB=65
PRIORITY_FILTER=70
PRIORITY_SITES=75
PRIORITY_GATHER=75
PRIORITY_UTIL=55

# 0.3 Used to find other files needed in repository
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# 0.4 Define and create job groups
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}"
JOB_GROUP_JOINT="/${USER}/compute-${COMPUTE_USER}/joint"

[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 500 ${JOB_GROUP}
[[ -z "$(bjgroup | grep $JOB_GROUP_JOINT)" ]] && bgadd -L 480 ${JOB_GROUP_JOINT}

## Begin Job submission
# Set up input
INPUT_VCF=$1
NAMEBASE=$(echo ${INPUT_VCF##*/} | cut -d. -f1)

# These 3 variables are used for each job submission to connect all the jobs for each sample consistent
JOBNAME="ngi-${USER}-${NAMEBASE}"
ENV_FILE="${SCRIPT_DIR}/../../baseEnvs/references_2_0.env"
LOGDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${NAMEBASE}
[[ -d $LOGDIR ]] || mkdir $LOGDIR

## 4. Gather data
# 4.1 Gather QCed vcfs
# This job gathers the filtered files in the provided directory.
# This job produces a gathered vcf file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-GATHER \
  -n 1 \
  -N \
  -o ${LOGDIR}/${NAMEBASE}.gather.%J.out \
  -R '{ select[mem>4GB] rusage[mem=4GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $(( PRIORITY_GATHER + 1 )) \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/gather_qced_vcfs.bash ${INPUT_VCF%/*}
