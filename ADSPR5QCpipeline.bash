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
PRIORITY_INTLIST=55
PRIORITY_MISS=60
PRIORITY_ANNAB=65
PRIORITY_FILTER=70
PRIORITY_SITES=75
PRIORITY_ANN=80
PRIORITY_GATHER=85
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
ENV_FILE="${SCRIPT_DIR}/baseEnvs/references_2_0.env"
LOGDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${NAMEBASE}
[[ -d $LOGDIR ]] || mkdir $LOGDIR

## 1. Prepare data
# 1.1 Make intlists
# This job creates intlists to use for scattering.
# This job produces intervals for pipeline use.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-INTLIST \
  -n 1 \
  -o ${LOGDIR}/${NAMEBASE}.intlist.%J.out \
  -Ne \
  -R '{ select[mem>4GB] rusage[mem=4GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_INTLIST \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/make_intlists.bash ${INPUT_VCF}

# 1.2 Set missing vcf
# This job sets genotypes with DP < 10 or GQ < 20 to missing.
# This job produces split files with a low quality genotypes set to missing vcf file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-SETMISSING[1-50] \
  -w "done(\"${JOBNAME}-INTLIST\")" \
  -n 1 \
  -o ${LOGDIR}/${NAMEBASE}.miss.%J.%I.out \
  -Ne \
  -R '{ select[mem>4GB] rusage[mem=4GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_MISS \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/set_vars_missing.bash ${INPUT_VCF}

## 2. Filter data
# 2.1 Annotate vcf with AB
# This job annotates the vcf with Allele Balance annotations.
# This job produces an annotated vcf file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-ANNAB[1-50] \
  -w "done(\"${JOBNAME}-SETMISSING[*]\")" \
  -n 1 \
  -o ${LOGDIR}/${NAMEBASE}.annAB.%J.%I.out \
  -Ne \
  -R '{ select[mem>20GB] rusage[mem=20GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_ANNAB \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/annotateAB_interval.bash ${INPUT_VCF%/*}

# 2.2 Filter vcf
# This job filters variants based on several metrics.
# This job produces a filtered vcf file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-FILTER[1-50] \
  -w "done(\"${JOBNAME}-ANNAB[*]\")" \
  -n 1 \
  -o ${LOGDIR}/${NAMEBASE}.filter.%J.%I.out \
  -Ne \
  -R '{ select[mem>4GB] rusage[mem=4GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_FILTER \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/GATKQC_filters.bash ${INPUT_VCF%/*}

# 3. Annotate data
# 3.1 Make sites-only vcf
# This job takes a filtered vcf and makes a sites-only vcf for annotation.
# This job produces a sites-only vcf file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-SITES[1-50] \
  -w "done(\"${JOBNAME}-FILTER[*]\")" \
  -n 1 \
  -o ${LOGDIR}/${NAMEBASE}.sites.%J.%I.out \
  -Ne \
  -R '{ select[mem>4GB] rusage[mem=4GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_SITES \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/make_sites_only_vcf.bash ${INPUT_VCF%/*}

# 3.2 Annotate vcf
# This job takes a sites-only vcf and annotates it with a variety of data.
# This job produces an annotated vcf file.
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_JOINT} \
  -J ${JOBNAME}-ANNOTATE[1-50] \
  -w "done(\"${JOBNAME}-SITES[*]\")" \
  -n 1 \
  -o ${LOGDIR}/${NAMEBASE}.ann.%J.%I.out \
  -Ne \
  -R '{ rusage[mem=40GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_ANN \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/annotateALL_interval.bash ${INPUT_VCF%/*}

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
  -w "done(\"${JOBNAME}-ANNOTATE\")" \
  -n 1 \
  -N \
  -o ${LOGDIR}/${NAMEBASE}.gather.%J.out \
  -R '{ select[mem>4GB] rusage[mem=4GB] }' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_GATHER \
  -a 'docker(mjohnsonngi/wxsjointqc:2.0)' \
  bash /scripts/gather_qced_vcfs.bash ${INPUT_VCF%/*}
