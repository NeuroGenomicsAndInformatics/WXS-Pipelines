#!/bin/bash
## Needed for Parabricks
export PATH="/opt/miniconda/bin:$PATH"

## Set up variables for specific users. These should be all that's needed to change user
export COMPUTE_USER=cruchagac
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
PRIORITY_ALIGN=60
PRIORITY_BQSR=65
PRIORITY_HC=70
PRIORITY_UTIL=80
PRIORITY_QC=50

# 0.3 Used to find other files needed in repository
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# 0.4 Define and create job groups
JOB_GROUP="/${USER}/compute-${COMPUTE_USER}"
JOB_GROUP_ALIGN="/${USER}/compute-${COMPUTE_USER}/align"
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"
[[ -z "$(bjgroup | grep $JOB_GROUP)" ]] && bgadd -L 300 ${JOB_GROUP}
[[ -z "$(bjgroup | grep $JOB_GROUP_ALIGN)" ]] && bgadd -L 20 ${JOB_GROUP_ALIGN}
[[ -z "$(bjgroup | grep $JOB_GROUP_QC)" ]] && bgadd -L 20 ${JOB_GROUP_QC}

## Begin Job submission loop
# This if checks if the submission was a workfile or a direct submission of a sample and adds them to an array
# Cram inputs require a workfile, but bams and paired end fastqs can be run with only the FULLSMID
if [[ -f $1 ]]; then FULLSMIDS=($(cat $1)); else FULLSMIDS=($@); fi
for FULLSMID in ${FULLSMIDS[@]}; do
# This script generates the .env environment file that all jobs use in the pipeline
# Each sample has a unique environment; however, the script adds the reference variables that are the same per version
bash ${SCRIPT_DIR}/makeSampleEnvWES.bash ${FULLSMID}
# These 3 variables are used for each job submission to connect all the jobs for each sample consistent
JOBNAME="ngi-${USER}-${FULLSMID}"
ENV_FILE="/scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/envs/$FULLSMID.env"
LOGDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}

## 1. Align, Sort, Intersect, and Mark Duplicates
# This job prepares data from all inputs to run as paired-end fastqs
# This job produces an aligned, sorted, duplicate-marked cram
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="${ENV_FILE}" \
bsub -g ${JOB_GROUP_ALIGN} \
  -J ${JOBNAME}-align \
  -n 8 \
  -o ${LOGDIR}/${FULLSMID}.align.%J.out \
  -R 'select[mem>80GB] rusage[mem=80GB/job] span[hosts=1]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -sp $PRIORITY_ALIGN \
  -a 'docker(mjohnsonngi/wxsalignerwes:2.0)' \
  bash /scripts/align.bash "$2"

## 2. BQSR
# This job produces a BQSR report
# The job itself runs BaseRecalibratorSpark on the cram created by the alignment job
# The generated BQSR report is used by haplotypecaller without having to do ApplyBQSR
LSF_DOCKER_VOLUMES="/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP} \
  -J ${JOBNAME}-bqsr \
  -w "done(\"${JOBNAME}-align\")" \
  -n 8 \
  -Ne \
  -sp $PRIORITY_BQSR \
  -o ${LOGDIR}/${FULLSMID}.bqsr.%J.out \
  -R 'select[mem>50GB] rusage[mem=50GB] span[hosts=1]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxsrecalibrator:2.0)' \
  bash /scripts/bqsrsparkwes.bash

## 3. Call Variants
# This job uses the cram from the alignment staged and the BQSR report to call variants
# This job produces a gvcf
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref \
$HOME:$HOME" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP} \
  -J ${JOBNAME}-hc \
  -w "done(\"${JOBNAME}-bqsr\")" \
  -n 2 \
  -Ne \
  -sp $PRIORITY_HC \
  -o ${LOGDIR}/${FULLSMID}.hc.%J.out \
  -R 'select[mem>50GB] rusage[mem=50GB]' \
  -G compute-${COMPUTE_USER} \
  -q general \
  -a 'docker(mjohnsonngi/wxshaplotypecallerwes:2.0)' \
  bash /scripts/cpuhcwes.bash

## 4. Stage out data
# This job just moves data from scratch to storage and cleans up
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
# 5.1 Coverage
# This job produces coverage reports for raw, high quality, and padded exome coverages
# This job uses the pipeline-generated cram while it's on Active storage
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-wgsmetrics \
    -w "done(\"${JOBNAME}-align\") && done(\"${JOBNAME}-stageout\")" \
    -n 2 \
    -Ne \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=25GB,tmp=2GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxscoverage:2.0)' \
    bash /scripts/gatkwgsmetrics_ex.bash ${FULLSMID}

# 5.2 FREEMIX
# This job uses VerifyBamID2 to create a report on possible contamination for the sample
# This job produces selfSM and Ancestry text files
# This job uses the pipeline-generated cram while it's on Active storage
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-freemix \
    -w "done(\"${JOBNAME}-align\") && done(\"${JOBNAME}-stageout\")" \
    -Ne \
    -n 2 \
    -sp $PRIORITY_QC \
    -R 'rusage[mem=20GB,tmp=2GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsfreemix:2.0)' \
    bash /scripts/vbid_exome.bash

## 5.3 Variant Calling Metrics
# This job produces a variant calling metrics report that includes Ti/Tv ratios and #s of SNPs and INDELS
# This job uses the pipeline-generated gvcf while it's on Active storage
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

## 5.4 Key Gene Annotations
# This job produces an annotation file using SnpEff
# The genes annotated include APP, PSEN1, PSEN2, GRN, TREM2, and MAPT
# This job uses the pipeline-generated gvcf while it's on Active storage
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

## 5.5 Stats File
# This job collects data from each of the previously generated reports into a single line
# This job produces a csv with a header line and a line of data from the QC reports
LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
/scratch1/fs1/${SCRATCH_USER}:/scratch1/fs1/${SCRATCH_USER} \
$HOME:$HOME \
$REF_DIR:/ref" \
LSF_DOCKER_ENV_FILE="$ENV_FILE" \
bsub -g ${JOB_GROUP_QC} \
    -J ${JOBNAME}-stageout \
    -w "ended(\"${JOBNAME}-wgsmetrics\") && ended(\"${JOBNAME}-vcfmetrics\") && ended(\"${JOBNAME}-freemix\") && ended(\"${JOBNAME}-snpeff\")" \
    -n 1 \
    -Ne \
    -sp $PRIORITY_UTIL \
    -o ${LOGDIR}/${FULLSMID}.stageout.%J.out \
    -R 'rusage[mem=4GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsstager:2.0)' \
    bash /scripts/statsupdate_exome.bash

done
