#!/bin/bash
ENVS_DIR="/scratch1/fs1/cruchagac/$USER/c1in/envs"
BASE_ENVS_DIR="./baseEnvs"
for FULLSMID in $(cat $1) ; do
ENV_FILE="$ENVS_DIR/${FULLSMID}.env"
echo -e "FULLSMID=${FULLSMID}" > $ENV_FILE
echo -e "INDIR=/input/${FULLSMID}" >> $ENV_FILE
mkdir /scratch1/fs1/cruchagac/$USER/c1out/$FULLSMID
echo -e "OUTDIR=/output/${FULLSMID}" >> $ENV_FILE
echo -e "LOGFILE=/output/${FULLSMID}/${FULLSMID}_runlog.txt" >> $ENV_FILE
echo -e "RUN_TYPE=paddedexome" >> $ENV_FILE
echo -e "RGBASES="$(basename -s .rgfile /scratch1/fs1/cruchagac/$USER/c1in/${FULLSMID}/*.rgfile)"" >> $ENVS_DIR/${FULLSMID}.env
cat ${BASE_ENVS_DIR}/pipelinebase.env >> $ENV_FILE
cat ${BASE_ENVS_DIR}/references.env >> $ENV_FILE
cp $ENV_FILE /scratch1/fs1/cruchagac/$USER/c1out/$FULLSMID/
done
