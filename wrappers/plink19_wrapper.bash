#!/bin/bash
REF_DIR="/scratch1/fs1/fernandezv/WXSref"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP="/${USER}/compute-fernandezv"

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
${REF_DIR}:/ref" \
bsub -g ${JOB_GROUP} \
    -J ngi-${USER}-plink19 \
    -n 8 \
    -N \
    -R 'rusage[mem=100GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1)' \
    $@
