#!/bin/bash
## Set up variables for specific users. These should be all that's needed to change user
export COMPUTE_USER=fernandezv
export SCRATCH_USER=cruchagac
export STORAGE_USER=cruchagac
export REF_DIR="/scratch1/fs1/cruchagac/WXSref"

## Touch the references so that compute1 doesn't remove them
find $REF_DIR -true -exec touch '{}' \;

## 0. Set up for job submission 
# 0.1 Make expected directories in case they are missing
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs

# 0.2 Priorities are set to handle bounded-buffer issues
PRIORITY_VQSR=60
PRIORITY_APPLY=65
PRIORITY_FILTER=70
PRIORITY_UTIL=55
PRIORITY_QC=50

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

# 1.4 Gather sites-only vcfs
# This job produces a gathered sites-only vcf file.
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
  -R '{ select[mem>20GB] rusage[mem=20GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_UTIL \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/gather_sites_only_vcfs.bash ${INPUT_VCF%/*}

## 2. VQSR
# 2.1 SNP VQSR
# This job runs the SNP VQSR recalibration task.
# This job produces a VQSR recal table and tranches file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-VQSR_SNP \
  -w "done(\"${JOBNAME}-GATHER\")" \
  -n 1 \
  -N \
  -o ${LOGDIR}/${NAMEBASE}.vqsrsnp.%J.out \
  -R '{ select[mem>50GB] rusage[mem=50GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_VQSR \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/VQSR_SNP.bash ${INPUT_VCF%.*.*}.ann.sites.gathered.vcf.gz

# 2.2 INDEL VQSR
# This job runs the INDEL VQSR recalibration task.
# This job produces a VQSR recal table and a tranches file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-VQSR_INDEL \
  -w "done(\"${JOBNAME}-GATHER\")" \
  -n 1 \
  -N \
  -o ${LOGDIR}/${NAMEBASE}.vqsrindel.%J.out \
  -R '{ select[mem>50GB] rusage[mem=50GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_VQSR \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/VQSR_INDEL.bash ${INPUT_VCF%.*.*}.ann.sites.gathered.vcf.gz

## 3. Apply VQSR
# 3.1 Apply SNP VQSR
# This job applies the SNP recal table to a vcf.
# This job is an array of jobs to scatter the work. The number of array jobs must match the number of scattered intervals.
LSF_DOCKER_VOLUMES="/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP} \
  -J ${JOBNAME}-Apply_VQSR_SNP[1-50] \
  -w "done(\"${JOBNAME}-VQSR_SNP\") && done(\"${JOBNAME}-VQSR_INDEL\")" \
  -n 1 \
  -Ne \
  -sp $PRIORITY_APPLY \
  -o ${LOGDIR}/${NAMEBASE}.applySNP.%J.%I.out \
  -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/ApplyVQSR.bash ${INPUT_VCF%/*} SNP ${INPUT_VCF%.*.*}.ann.sites.gathered.SNP_recalibrate.recal

# 3.2 Apply INDEL VQSR
# This job applies the INDEL recal table to a vcf.
# This job is an array of jobs to scatter the work. The number of array jobs must match the number of scattered intervals.
LSF_DOCKER_VOLUMES="/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP} \
  -J ${JOBNAME}-Apply_VQSR_INDEL[1-50] \
  -w "done(\"${JOBNAME}-Apply_VQSR_SNP[*]\")" \
  -n 1 \
  -Ne \
  -sp $PRIORITY_APPLY \
  -o ${LOGDIR}/${NAMEBASE}.applyINDEL.%J.%I.out \
  -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/ApplyVQSR.bash ${INPUT_VCF%/*} INDEL ${INPUT_VCF%.*.*}.ann.sites.gathered.INDEL_recalibrate.recal

## 4. GATK QC Filtering
# This job runs the hard filtering on an interval vcf.
# This job produces a filtered vcf.
LSF_DOCKER_VOLUMES="/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP} \
  -J ${JOBNAME}-HARDFILT[1-50] \
  -w "done(\"${JOBNAME}-Apply_VQSR_INDEL[*]\")" \
  -n 1 \
  -Ne \
  -sp $PRIORITY_FILTER \
  -o ${LOGDIR}/${NAMEBASE}.hardfilt.%J.%I.out \
  -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/HardFilter.bash ${INPUT_VCF%.*.*}.ann.gathered.SNP.INDEL.vcf.gz
