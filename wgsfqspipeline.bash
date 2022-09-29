#!/bin/bash
## Needed for Parabricks
export PATH="/opt/miniconda/bin:$PATH"

## Set up variables
REF_DIR="/scratch1/fs1/cruchagac/WXSref"
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
JOB_GROUP_GPU="/${USER}-gpu/compute-fernandezv"
JOB_GROUP_ALIGN="/${USER}/compute-fernandezv/align"
[[ -z "$(bjgroup | grep $JOB_GROUP_C)" ]] && bgadd -L 10 ${JOB_GROUP_C}
[[ -z "$(bjgroup | grep $JOB_GROUP_F)" ]] && bgadd -L 300 ${JOB_GROUP_F}
[[ -z "$(bjgroup | grep $JOB_GROUP_GPU)" ]] && bgadd -L 15 ${JOB_GROUP_GPU}
[[ -z "$(bjgroup | grep $JOB_GROUP_ALIGN)" ]] && bgadd -L 15 ${JOB_GROUP_ALIGN}

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
bsub -g ${JOB_GROUP_GPU} \
  -J ${JOBNAME}-aligngpu \
  -n "1,16" \
  -Ne \
  -o ${LOGDIR}/${FULLSMID}.fq2bam.%J.out \
  -R '{ 16*{ select[gpuhost && mem>180GB] rusage[mem=180GB/job, ngpus_physical=1:gmodel=NVIDIAA100_SXM4_40GB] span[hosts=1] } } || { 16*{ select[gpuhost && mem>180GB] rusage[mem=180GB/job, ngpus_physical=1:gmodel=TeslaV100_SXM2_32GB] span[hosts=1] } }@2 || { 1*{ select[!gpuhost] } }@10' \
  -G compute-fernandezv \
  -q general \
  -sp $PRIORITY_ALIGN \
  -a 'docker(mjohnsonngi/wxsaligner:2.0)' \
  bash /scripts/align.bash

  LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
  /scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
  /scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
  /scratch1/fs1/ris/application/parabricks-license:/opt/parabricks \
  ${REF_DIR}:/ref \
  $HOME:$HOME" \
  LSF_DOCKER_NETWORK=host \
  LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
  LSF_DOCKER_ENTRYPOINT=/bin/sh \
  LSF_DOCKER_ENV_FILE="$ENV_FILE" \
  bsub -g ${JOB_GROUP_ALIGN} \
    -J ${JOBNAME}-aligncpu \
    -w "exit(\"${JOBNAME}-aligngpu\",66)" \
    -n 1 \
    -Ne \
    -o ${LOGDIR}/${FULLSMID}.fq2bam.%J.out \
    -R 'rusage[mem=10GB]' \
    -G compute-fernandezv \
    -q general \
    -sp $PRIORITY_ALIGN \
    -a 'docker(mjohnsonngi/wxsaligner:2.0)' \
    bash /scripts/stageinfqsalign3.bash

## 2. BQSR
LSF_DOCKER_VOLUMES="/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
  -J ${JOBNAME}-bqsr \
  -w "done(\"${JOBNAME}-aligngpu\") || done(\"${JOBNAME}-aligncpu\")" \
  -n 8 \
  -Ne \
  -sp $PRIORITY_BQSR \
  -o ${LOGDIR}/${FULLSMID}.bqsr.%J.out \
  -R 'select[mem>120GB] rusage[mem=120GB] span[hosts=1]' \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxsrecalibrator:2.0)' \
  bash /scripts/bqsrspark.bash

## 3. Call Variants
LSF_DOCKER_VOLUMES="/scratch1/fs1/fernandezv:/scratch1/fs1/fernandezv \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
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
  -n 16 \
  -Ne \
  -sp $PRIORITY_HC \
  -o ${LOGDIR}/${FULLSMID}.hc.%J.out \
  -R 'select[gpuhost && mem>140GB] rusage[mem=140GB] span[hosts=1]' \
  -gpu "num=1:j_exclusive=yes" \
  -G compute-fernandezv \
  -q general \
  -a 'docker(mjohnsonngi/wxshaplotypecaller:2.0)' \
  bash /scripts/gpuhc.bash

## 4. Stage out data
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
$HOME:$HOME" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-stageout \
    -w "exit(\"${JOBNAME}-bqsr\") || ended(\"${JOBNAME}-hc\")" \
    -n 4 \
    -sp $PRIORITY_UTIL \
    -o ${LOGDIR}/${FULLSMID}.stageout.%J.out \
    -R 'rusage[mem=4GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxspipeline:1.1)' \
    bash $SCRIPT_DIR/stageoutcram.bash $FULLSMID

## 5. QC
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-wgsmetrics \
    -w "(done(\"${JOBNAME}-aligngpu\") || done(\"${JOBNAME}-aligncpu\")) && done(\"${JOBNAME}-stageout\")" \
    -n 2 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=80GB,tmp=2GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxscoverage:2.0)' \
    bash /scripts/get_both_wgsmetrics.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-freemix \
    -w "(done(\"${JOBNAME}-aligngpu\") || done(\"${JOBNAME}-aligncpu\")) && done(\"${JOBNAME}-stageout\")" \
    -n 4 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=80GB,tmp=2GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxsfreemix:2.0)' \
    bash /scripts/vbid.bash ${FINAL_OUTDIR}/${CRAM}

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-vcfmetrics \
    -w "done(\"${JOBNAME}-hc\") && done(\"${JOBNAME}-stageout\")" \
    -n 8 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=40GB,tmp=2GB]' \
    -G compute-cruchagac \
    -q general \
    -a 'docker(mjohnsonngi/wxsvariantmetrics:2.0)' \
    bash /scripts/getvcfmetrics.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
${REF_DIR}:/ref" \
LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_F} \
    -J ${JOBNAME}-snpeff \
    -w "done(\"${JOBNAME}-hc\") && done(\"${JOBNAME}-stageout\")" \
    -n 16 \
    -sp $PRIORITY_HC \
    -o ${LOGDIR}/${FULLSMID}.snpeff.%J.out \
    -R 'rusage[mem=120GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxskeygeneannotator:2.0)' \
  	bash /scripts/keygene_annotate.bash
done
