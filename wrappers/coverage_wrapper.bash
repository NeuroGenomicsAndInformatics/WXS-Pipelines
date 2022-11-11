#!/bin/bash
REF_DIR="/scratch1/fs1/fernandezv/WXSref"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-fernandezv/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-wgsmetrics \
    -n 2 \
    -Ne \
    -R 'rusage[mem=25GB,tmp=10GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxscoverage:2.0)' \
    bash /scripts/get_all_wgsmetrics.bash $1
