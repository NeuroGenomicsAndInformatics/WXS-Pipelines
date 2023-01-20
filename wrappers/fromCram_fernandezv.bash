#!/bin/bash
## Needed for Parabricks
export PATH="/opt/miniconda/bin:$PATH"
REF_DIR="/scratch1/fs1/fernandezv/WXSref"

[ ! -d /scratch1/fs1/fernandezv/${USER} ] && mkdir /scratch1/fs1/fernandezv/${USER}
[ ! -d /scratch1/fs1/fernandezv/${USER}/c1in ] && mkdir /scratch1/fs1/fernandezv/${USER}/c1in
[ ! -d /scratch1/fs1/fernandezv/${USER}/c1out ] && mkdir /scratch1/fs1/fernandezv/${USER}/c1out
[ ! -d /storage1/fs1/cruchagac/Active/${USER}/c1out ] && mkdir //storage1/fs1/cruchagac/Active/${USER}/c1out
[ ! -d /scratch1/fs1/fernandezv/${USER}/c1out/logs ] && mkdir /scratch1/fs1/fernandezv/${USER}/c1out/logs
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
bash $HOME/WXS-Pipelines/makeSampleEnv_vf.bash ${FULLSMID}
JOBNAME="ngi-${USER}-${FULLSMID}"
ENV_FILE=/scratch1/fs1/fernandezv/${USER}/c1in/envs/${FULLSMID}.env
LOGDIR=/scratch1/fs1/fernandezv/${USER}/c1out/logs/${FULLSMID}

## 2. BQSR
LSF_DOCKER_VOLUMES="/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
  -J ${JOBNAME}-bqsr \
  -n 8 \
  -Ne \
  -sp $PRIORITY_BQSR \
  -o ${LOGDIR}/${FULLSMID}.bqsr.%J.out \
  -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]' \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxsrecalibrator:2.0)' \
  bash $SCRIPT_DIR/fixstage.bash\; \
  bash /scripts/bqsrspark.bash

## 3. Call Variants
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
/scratch1/fs1/ris/application/parabricks:/opt/parabricks \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/sh \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP_GPU} \
  -J ${JOBNAME}-hc \
  -w "done(\"${JOBNAME}-bqsr\")" \
  -n 8 \
  -Ne \
  -sp $PRIORITY_HC \
  -o ${LOGDIR}/${FULLSMID}.hc.%J.out \
  -R 'select[gpuhost && mem>180GB] rusage[mem=180GB] span[hosts=1]' \
  -gpu "num=1:gmem=16GB:j_exclusive=yes" \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxshaplotypecaller:2.0)' \
  bash /scripts/gpuhc.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE=$ENV_FILE \
bsub -g ${JOB_GROUP_F} \
  -J ${JOBNAME}-stageout \
  -w "ended(\"${JOBNAME}-hc\")" \
  -n1 \
  -sp 90 \
  -R 'rusage[mem=4GB]' \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxsstager:2.0)' \
  bash /scripts/stageout.bash
done
