#!/bin/bash
## Needed for Parabricks
export PATH="/opt/miniconda/bin:$PATH"

## Set up variables and update touch references
STORAGE_REF_DIR="/storage1/fs1/cruchagac/Active/matthew.j/REF/WXSref"
REF_DIR="/scratch1/fs1/cruchagac/WXSref"
find $REF_DIR -true -exec touch '{}' \;

[ ! -d /scratch1/fs1/cruchagac/${USER} ] && mkdir /scratch1/fs1/cruchagac/${USER}
[ ! -d /scratch1/fs1/cruchagac/${USER}/c1in ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1in
[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1out
[ ! -d /storage1/fs1/cruchagac/Active/${USER}/c1out ] && mkdir //storage1/fs1/cruchagac/Active/${USER}/c1out
[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/logs ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1out/logs
PRIORITY_ALIGN=60
PRIORITY_BQSR=65
PRIORITY_HC=70
PRIORITY_UTIL=80
PRIORITY_QC=50

## Set up directories and job submission variables
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_C="/${USER}/compute-cruchagac"
JOB_GROUP_F="/${USER}/compute-fernandezv"
JOB_GROUP_GPU="/${USER}/compute-fernandezv/gpu"
JOB_GROUP_ALIGN="/${USER}/compute-fernandezv/align"
[[ -z "$(bjgroup | grep $JOB_GROUP_C)" ]] && bgadd -L 10 ${JOB_GROUP_C}
[[ -z "$(bjgroup | grep $JOB_GROUP_F)" ]] && bgadd -L 300 ${JOB_GROUP_F}
[[ -z "$(bjgroup | grep $JOB_GROUP_GPU)" ]] && bgadd -L 10 ${JOB_GROUP_GPU}
[[ -z "$(bjgroup | grep $JOB_GROUP_ALIGN)" ]] && bgadd -L 20 ${JOB_GROUP_ALIGN}

if [[ -f $1 ]]; then FULLSMIDS=($(cat $1)); else FULLSMIDS=($@); fi
for FULLSMID in ${FULLSMIDS[@]}; do
bash ${SCRIPT_DIR}/makeSampleEnv.bash ${FULLSMID}
JOBNAME="ngi-${USER}-${FULLSMID}"
ENV_FILE="/scratch1/fs1/cruchagac/${USER}/c1in/envs/$FULLSMID.env"
LOGDIR=/scratch1/fs1/cruchagac/${USER}/c1out/logs/${FULLSMID}

## 1. Align
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
/scratch1/fs1/ris/application/parabricks-license:/opt/parabricks \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/bash \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_ALIGN} \
  -J ${JOBNAME}-aligngpu \
  -n12 \
  -o ${LOGDIR}/${FULLSMID}.fq2bam.%J.out \
  -R 'select[gpuhost && mem>120GB] rusage[mem=120GB/job, ngpus_physical=2:gmem=24GB] span[hosts=1]' \
  -G compute-fernandezv \
  -q general \
  -sp $PRIORITY_ALIGN \
  -a 'docker(mjohnsonngi/wxsaligner:2.0)' \
  bash /scripts/gpualign2.bash

done
