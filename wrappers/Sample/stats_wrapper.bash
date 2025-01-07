#!/bin/bash
# This wrapper generates a stats.csv file from the various QC reports
# The argument for this wrapper is a path to the sample's output directory on Active storage
FINAL_OUTDIR=$1

STORAGE_USER=cruchagac
COMPUTE_USER=fernandezv
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-stats \
    -n 1 \
    -Ne \
    -R 'rusage[mem=25GB,tmp=10GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsstager:2.0)' \
    bash /scripts/statsupdate.bash $FINAL_OUTDIR
