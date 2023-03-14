#!/bin/bash
## This wrapper script takes the data from a failed pipeline run and tries to generate a gvcf
# This wrapper assumes the cram and BQSR report are in the output directory on Active storage
# The argument for this wrapper script is a FULLSMID or workfile
# This script is essentially a subset of the job submissions from the pipeline script

## Needed for Parabricks
export PATH="/opt/miniconda/bin:$PATH"

REF_DIR="/scratch1/fs1/fernandezv/WXSref"
export COMPUTE_USER=fernandezv
export SCRATCH_USER=cruchagac
export STORAGE_USER=cruchagac

## Set up directories and job submission variables
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs
PRIORITY_ALIGN=60
PRIORITY_BQSR=65
PRIORITY_HC=70
PRIORITY_UTIL=80
PRIORITY_QC=50

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}"
JOB_GROUP_GPU="/${USER}/compute-${COMPUTE_USER}/gpu"
JOB_GROUP_ALIGN="/${USER}/compute-${COMPUTE_USER}/align"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 300 ${JOB_GROUP}
[[ -z "$(bjgroup | grep $JOB_GROUP_GPU)" ]] && bgadd -L 10 ${JOB_GROUP_GPU}
[[ -z "$(bjgroup | grep $JOB_GROUP_ALIGN)" ]] && bgadd -L 20 ${JOB_GROUP_ALIGN}

if [[ -f $1 ]]; then FULLSMIDS=($(cat $1)); else FULLSMIDS=($@); fi
for FULLSMID in ${FULLSMIDS[@]}; do
bash $SCRIPT_DIR/../makeSampleEnv.bash ${FULLSMID}
JOBNAME="ngi-${USER}-${FULLSMID}"
ENV_FILE=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/envs/${FULLSMID}.env
LOGDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}

## 3. Call Variants
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
/scratch1/fs1/ris/application/parabricks:/opt/parabricks \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/sh \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP_GPU} \
  -J ${JOBNAME}-hc \
  -n 8 \
  -Ne \
  -sp $PRIORITY_HC \
  -o ${LOGDIR}/${FULLSMID}.hc.%J.out \
  -R 'select[gpuhost && mem>180GB] rusage[mem=180GB] span[hosts=1]' \
  -gpu "num=1:gmem=16GB:j_exclusive=yes" \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxshaplotypecaller:2.0)' \
  bash $SCRIPT_DIR/fixstage.bash\; \
  bash /scripts/gpuhc.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP} \
  -J ${JOBNAME}-stageout \
  -w "ended(\"${JOBNAME}-hc\")" \
  -n1 \
  -sp 90 \
  -R 'rusage[mem=4GB]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxsstager:2.0)' \
  bash /scripts/stageout.bash
done
