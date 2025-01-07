#!/bin/bash
# This wrapper script generates the annotation file with SnpEff
# The genes annotated include PSEN1, PSEN2, APP, GRN, TREM2, and MAPT
# The argument for this wrapper script is the full path to a gvcf
FULL_GVCF=$1

STORAGE_USER=cruchagac
COMPUTE_USER=fernandezv
REF_DIR="/scratch1/fs1/cruchagac/WXSref"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
JOB_GROUP_QC="/${USER}/compute-${COMPUTE_USER}/qc"

LSF_DOCKER_VOLUMES="/storage1/fs1/${STORAGE_USER}/Active:/storage1/fs1/${STORAGE_USER}/Active \
${REF_DIR}:/ref" \
LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
LSF_DOCKER_ENV_FILE=$SCRIPT_DIR/../../baseEnvs/references_2_0.env \
bsub -g ${JOB_GROUP_QC} \
    -J ngi-${USER}-snpeff \
    -Ne \
    -n 2 \
    -R 'rusage[mem=25GB]' \
    -G compute-${COMPUTE_USER} \
    -q general \
    -a 'docker(mjohnsonngi/wxskeygeneannotator:2.0)' \
  	bash /scripts/allgene_annotate.bash $FULL_GVCF
