#!/bin/bash
## Needed for Parabricks
export PATH="/opt/miniconda/bin:$PATH"

## Set up variables and update touch references
REF_DIR="/scratch1/fs1/fernandezv/WXSref"
find $REF_DIR -true -exec touch '{}' \;

[ ! -d /scratch1/fs1/fernandezv/${USER} ] && mkdir /scratch1/fs1/fernandezv/${USER}
[ ! -d /scratch1/fs1/fernandezv/${USER}/c1in ] && mkdir /scratch1/fs1/fernandezv/${USER}/c1in
[ ! -d /scratch1/fs1/fernandezv/${USER}/c1out ] && mkdir /scratch1/fs1/fernandezv/${USER}/c1out
[ ! -d /storage1/fs1/cruchagac/Active/${USER}/c1out ] && mkdir /storage1/fs1/cruchagac/Active/${USER}/c1out
[ ! -d /scratch1/fs1/fernandezv/${USER}/c1out/logs ] && mkdir /scratch1/fs1/fernandezv/${USER}/c1out/logs
PRIORITY_ALIGN=60
PRIORITY_BQSR=65
PRIORITY_HC=70
PRIORITY_UTIL=80
PRIORITY_QC=50

## Set up directories and job submission variables
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_F="/${USER}/compute-fernandezv"
JOB_GROUP_GPU="/${USER}/compute-fernandezv/gpu"
JOB_GROUP_ALIGN="/${USER}/compute-fernandezv/align"
[[ -z "$(bjgroup | grep $JOB_GROUP_C)" ]] && bgadd -L 10 ${JOB_GROUP_C}
[[ -z "$(bjgroup | grep $JOB_GROUP_F)" ]] && bgadd -L 300 ${JOB_GROUP_F}
[[ -z "$(bjgroup | grep $JOB_GROUP_GPU)" ]] && bgadd -L 10 ${JOB_GROUP_GPU}
[[ -z "$(bjgroup | grep $JOB_GROUP_ALIGN)" ]] && bgadd -L 20 ${JOB_GROUP_ALIGN}

if [[ -f $1 ]]; then FULLSMIDS=($(cat $1)); else FULLSMIDS=($@); fi
for FULLSMID in ${FULLSMIDS[@]}; do
bash ${SCRIPT_DIR}/makeSampleEnv_vf.bash ${FULLSMID}
JOBNAME="ngi-${USER}-${FULLSMID}"
ENV_FILE="/scratch1/fs1/fernandezv/${USER}/c1in/envs/$FULLSMID.env"
LOGDIR=/scratch1/fs1/fernandezv/${USER}/c1out/logs/${FULLSMID}

## 1. Align
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
/scratch1/fs1/ris/application/parabricks-license:/opt/parabricks \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/bash \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_ALIGN} \
  -J ${JOBNAME}-align \
  -n8 \
  -o ${LOGDIR}/${FULLSMID}.fq2bam.%J.out \
  -R '{ select[gpuhost && mem>180GB] rusage[ngpus_physical=1:gmem=12GB, mem=180GB/job] span[hosts=1] } || { select[!gpuhost] rusage[mem=180GB/job] }@10' \
  -G compute-fernandezv \
  -q general \
  -sp $PRIORITY_ALIGN \
  -a 'docker(mjohnsonngi/wxsaligner:2.0)' \
  bash /scripts/align.bash

## 2. BQSR
LSF_DOCKER_VOLUMES="/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
  -J ${JOBNAME}-bqsr \
  -w "done(\"${JOBNAME}-align\")" \
  -n 8 \
  -Ne \
  -sp $PRIORITY_BQSR \
  -o ${LOGDIR}/${FULLSMID}.bqsr.%J.out \
  -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]' \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxsrecalibrator:2.0)' \
  bash /scripts/bqsrspark.bash

## 3. Call Variants
LSF_DOCKER_VOLUMES="/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
/scratch1/fs1/ris/application/parabricks:/opt/parabricks \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/sh \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_GPU} \
  -J ${JOBNAME}-hc \
  -w "done(\"${JOBNAME}-bqsr\")" \
  -n 8 \
  -Ne \
  -sp $PRIORITY_HC \
  -o ${LOGDIR}/${FULLSMID}.hc.%J.out \
  -R 'select[gpuhost && mem>140GB] rusage[mem=140GB] span[hosts=1]' \
  -gpu "num=1:gmem=12GB:j_exclusive=yes" \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxshaplotypecaller:2.0)' \
  bash /scripts/gpuhc.bash

## 4. Stage out data
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-stageout \
    -w "exit(\"${JOBNAME}-bqsr\") || ended(\"${JOBNAME}-hc\")" \
    -n 1 \
    -sp $PRIORITY_UTIL \
    -o ${LOGDIR}/${FULLSMID}.stageout.%J.out \
    -R 'rusage[mem=4GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxsstager:2.0)' \
    bash /scripts/stageout.bash

## 5. QC
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-wgsmetrics \
    -w "done(\"${JOBNAME}-align\") && done(\"${JOBNAME}-stageout\")" \
    -n 2 \
    -Ne \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=25GB,tmp=2GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxscoverage:2.0)' \
    bash /scripts/get_all_wgsmetrics.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-freemix \
    -w "done(\"${JOBNAME}-align\") && done(\"${JOBNAME}-stageout\")" \
    -Ne \
    -n 4 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=80GB,tmp=2GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxsfreemix:2.0)' \
    bash /scripts/vbid.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-vcfmetrics \
    -w "done(\"${JOBNAME}-hc\") && done(\"${JOBNAME}-stageout\")" \
    -Ne \
    -n 8 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=40GB,tmp=2GB]' \
    -G compute-cruchagac \
    -q general \
    -a 'docker(mjohnsonngi/wxsvariantmetrics:2.0)' \
    bash /scripts/gatkvcfmetrics.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
${REF_DIR}:/ref" \
LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-snpeff \
    -w "done(\"${JOBNAME}-hc\") && done(\"${JOBNAME}-stageout\")" \
    -Ne \
    -n 4 \
    -sp $PRIORITY_HC \
    -o ${LOGDIR}/${FULLSMID}.snpeff.%J.out \
    -R 'rusage[mem=120GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxskeygeneannotator:2.0)' \
  	bash /scripts/keygene_annotate.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-stageout \
    -w "ended(\"${JOBNAME}-wgsmetrics\") && ended(\"${JOBNAME}-vcfmetrics\") && ended(\"${JOBNAME}-freemix\") && ended(\"${JOBNAME}-snpeff\")" \
    -n 1 \
    -N \
    -sp $PRIORITY_UTIL \
    -o ${LOGDIR}/${FULLSMID}.stageout.%J.out \
    -R 'rusage[mem=4GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxsstager:2.0)' \
    bash /scripts/statsupdate.bash

done