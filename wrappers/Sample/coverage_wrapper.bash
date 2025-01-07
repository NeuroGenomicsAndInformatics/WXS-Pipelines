#!/bin/bash
# This wrapper generates raw, high quality, and padded exome coverage reports
# The argument is the full path to the cram to be analyzed
FULL_CRAM=$1

STORAGE_USER=cruchagac
COMPUTE_USER=fernandezv
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-wgsmetrics \
    -n 2 \
    -Ne \
    -R 'rusage[mem=25GB,tmp=10GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxscoverage:2.0)' \
    bash /scripts/get_all_wgsmetrics.bash $FULL_CRAM
