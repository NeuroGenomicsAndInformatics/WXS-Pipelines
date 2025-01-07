#!/bin/bash
# This wrapper runs variant calling metrics that gives data like ti/tv ratio
# The argument is the full path to a gvcf for the sample
FULL_GVCF=$1

STORAGE_USER=cruchagac
COMPUTE_USER=fernandezv
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-vcfmetrics \
    -Ne \
    -n 4 \
    -R 'rusage[mem=10GB,tmp=2GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsvariantmetrics:2.0)' \
    bash /scripts/gatkvcfmetrics.bash $FULL_GVCF
