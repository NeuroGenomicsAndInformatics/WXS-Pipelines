#!/bin/bash
CRAM=$1
OUTDIR=$2

REF_DIR="/scratch1/fs1/fernandezv/WXSref"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-fernandezv/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac:/storage1/fs1/cruchagac \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-archive \
    -n 4 \
    -Ne \
    -R 'rusage[mem=25GB,tmp=10GB]' \
    -G compute-fernandezv \
    -q general \
    -a 'docker(mjohnsonngi/wxsarchiver:2.0)' \
    bash /scripts/archive_cram.bash $CRAM $OUTDIR
