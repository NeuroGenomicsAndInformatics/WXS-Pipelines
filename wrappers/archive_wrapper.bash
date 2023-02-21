#!/bin/bash
# This wrapper script generates an archived cram using samtools
# The first argument is the full path to a cram 
# The second argument is a path to the directory where the archived cram will go
# The second argument shouldn't end in /
CRAM=$1
OUTDIR=$2

STORAGE_USER=cruchagac
COMPUTE_USER=fernandezv
REF_DIR="/scratch1/fs1/fernandezv/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}:/storage1/fs1/${STORAGE_USER} \
${REF_DIR}:/ref" \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-archive \
    -n 4 \
    -Ne \
    -R 'rusage[mem=25GB,tmp=10GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxsarchiver:2.0)' \
    bash /scripts/archive_cram.bash $CRAM $OUTDIR
