#!/bin/bash
ENVS_DIR="/scratch1/fs1/cruchagac/${USER}/c1in/envs"
BASE_ENVS_DIR="./baseEnvs"
LINE="$1"
FULLSMID="$(echo $LINE | cut -d ' ' -f1)"
export CRAM="$(echo $LINE | cut -d ' ' -f2)"
STAGE1_WORKFILE="~/work${CRAM##/*}.txt"
ENV_FILE="$ENVS_DIR/${CRAM}.env"
echo -e "FULLSMID=${FULLSMID}" > $ENV_FILE
echo -e "STAGE_INDIR=/storage1/fs1/cruchagac/Active/${USER}/c1in" >> $ENV_FILE
echo -e "INDIR=/input/${FULLSMID}" >> $ENV_FILE
[ ! -d /scratch1/fs1/cruchagac/$USER/c1out/${FULLSMID} ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1out/$FULLSMID
echo -e "OUTDIR=/output/${FULLSMID}" >> $ENV_FILE
echo -e "STAGE1_WORKFILE=${STAGE1_WORKFILE}" >> $ENV_FILE
cat ${BASE_ENVS_DIR}/pipelinebase.env >> $ENV_FILE
cat ${BASE_ENVS_DIR}/references.env >> $ENV_FILE
rsync $ENV_FILE /scratch1/fs1/cruchagac/${USER}/c1out/$FULLSMID/
