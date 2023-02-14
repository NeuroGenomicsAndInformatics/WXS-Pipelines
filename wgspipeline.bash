#!/bin/bash
## Needed for Parabricks
export PATH="/opt/miniconda/bin:$PATH"

## Set up variables and update touch references
export COMPUTE_USER=fernandezv
export SCRATCH_USER=cruchagac
export STORAGE_USER=cruchagac
STORAGE_REF_DIR="/storage1/fs1/cruchagac/Active/matthew.j/REF/WXSref"
export REF_DIR="/scratch1/fs1/fernandezv/WXSref"
find $REF_DIR -true -exec touch '{}' \;

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

## Set up directories and job submission variables
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}"
JOB_GROUP_GPU="/${USER}/compute-${COMPUTE_USER}/gpu"
JOB_GROUP_ALIGN="/${USER}/compute-${COMPUTE_USER}/align"
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 300 ${JOB_GROUP_F}
[[ -z "$(bjgroup | grep $JOB_GROUP_GPU)" ]] && bgadd -L 10 ${JOB_GROUP_GPU}
[[ -z "$(bjgroup | grep $JOB_GROUP_ALIGN)" ]] && bgadd -L 20 ${JOB_GROUP_ALIGN}
[[ -z "$(bjgroup | grep $JOB_GROUP_QC)" ]] && bgadd -L 20 ${JOB_GROUP_QC}

if [[ -f $1 ]]; then FULLSMIDS=($(cat $1)); else FULLSMIDS=($@); fi
for FULLSMID in ${FULLSMIDS[@]}; do
bash ${SCRIPT_DIR}/makeSampleEnv.bash ${FULLSMID}
JOBNAME="ngi-${USER}-${FULLSMID}"
ENV_FILE="/scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/envs/$FULLSMID.env"
LOGDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}

## 1. Align
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
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
  -R '{ select[gpuhost && mem>220GB] rusage[ngpus_physical=1:gmem=12GB, mem=220GB/job] span[hosts=1] } || { select[!gpuhost] rusage[mem=220GB/job] }@10' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_ALIGN \
  -a 'docker(mjohnsonngi/wxsaligner:2.0)' \
  bash /scripts/align.bash "$2"

## Fallback if GPU fails
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
/scratch1/fs1/ris/application/parabricks-license:/opt/parabricks \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/bash \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_ALIGN} \
  -J ${JOBNAME}-align2 \
  -w "exit(\"${JOBNAME}-align\")" \
  -n8 \
  -o ${LOGDIR}/${FULLSMID}.fq2bam.%J.out \
  -R 'select[mem>180GB] rusage[mem=180GB/job] span[hosts=1]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_ALIGN \
  -a 'docker(mjohnsonngi/wxsaligner:2.0)' \
  bash /scripts/align.bash "$2"

## 2. BQSR
LSF_DOCKER_VOLUMES="/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP} \
  -J ${JOBNAME}-bqsr \
  -w "done(\"${JOBNAME}-align\") || done(\"${JOBNAME}-align2\")" \
  -n 8 \
  -Ne \
  -sp $PRIORITY_BQSR \
  -o ${LOGDIR}/${FULLSMID}.bqsr.%J.out \
  -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxsrecalibrator:2.0)' \
  bash /scripts/bqsrspark.bash

## 3. Call Variants
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
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
  -R 'select[gpuhost && mem>180GB] rusage[mem=180GB] span[hosts=1]' \
  -gpu "num=1:gmem=12GB:j_exclusive=yes" \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxshaplotypecaller:2.0)' \
  bash /scripts/gpuhc.bash

## 4. Stage out data
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP} \
    -J ${JOBNAME}-stageout \
    -w "exit(\"${JOBNAME}-bqsr\") || ended(\"${JOBNAME}-hc\")" \
    -n 1 \
    -sp $PRIORITY_UTIL \
    -o ${LOGDIR}/${FULLSMID}.stageout.%J.out \
    -R 'rusage[mem=4GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsstager:2.0)' \
    bash /scripts/stageout.bash\; sleep 100

## 5. QC
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-wgsmetrics \
    -w "(done(\"${JOBNAME}-align\") || done(\"${JOBNAME}-align2\")) && done(\"${JOBNAME}-stageout\")" \
    -n 2 \
    -Ne \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=25GB,tmp=2GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxscoverage:2.0)' \
    bash /scripts/get_all_wgsmetrics.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-freemix \
    -w "(done(\"${JOBNAME}-align\") || done(\"${JOBNAME}-align2\")) && done(\"${JOBNAME}-stageout\")" \
    -Ne \
    -n 2 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=20GB,tmp=2GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsfreemix:2.0)' \
    bash /scripts/vbid.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-vcfmetrics \
    -w "done(\"${JOBNAME}-hc\") && done(\"${JOBNAME}-stageout\")" \
    -Ne \
    -n 4 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=10GB,tmp=2GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsvariantmetrics:2.0)' \
    bash /scripts/gatkvcfmetrics.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref" \
LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-snpeff \
    -w "done(\"${JOBNAME}-hc\") && done(\"${JOBNAME}-stageout\")" \
    -Ne \
    -n 2 \
    -sp $PRIORITY_QC \
    -o ${LOGDIR}/${FULLSMID}.snpeff.%J.out \
    -R 'rusage[mem=25GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxskeygeneannotator:2.0)' \
  	bash /scripts/keygene_annotate.bash

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-stageout \
    -w "ended(\"${JOBNAME}-wgsmetrics\") && ended(\"${JOBNAME}-vcfmetrics\") && ended(\"${JOBNAME}-freemix\") && ended(\"${JOBNAME}-snpeff\")" \
    -n 1 \
    -N \
    -sp $PRIORITY_UTIL \
    -o ${LOGDIR}/${FULLSMID}.stageout.%J.out \
    -R 'rusage[mem=4GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsstager:2.0)' \
    bash /scripts/statsupdate.bash

done
