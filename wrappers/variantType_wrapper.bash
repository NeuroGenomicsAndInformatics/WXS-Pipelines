#!/bin/bash
# This wrapper generates a log file containing variants split by type
# The argument for this wrapper is a vcf file
VCF=$1

SCRATCH_USER=fernandezv
STORAGE_USER=cruchagac
COMPUTE_USER=fernandezv
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-varianttype \
    -n 1 \
    -Ne \
    -R 'rusage[mem=4GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsjointasqc:2.0)' \
    bash /scripts/variantTypes.bash $VCF
