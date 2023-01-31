#!/bin/bash
REF_DIR="/scratch1/fs1/fernandezv/WXSref"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP="/${USER}/compute-fernandezv"

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
${REF_DIR}:/ref" \
bsub -g ${JOB_GROUP} \
    -J ngi-${USER}-sra \
    -n 1 \
    -Ne \
    -R 'rusage[mem=10GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(pegi3s/sratoolkit:latest)' \
    $@
