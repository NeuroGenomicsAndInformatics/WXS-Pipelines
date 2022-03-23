#!/bin/bash
ENVS_DIR="/scratch1/fs1/cruchagac/${USER}/c1in/envs"
BASE_ENVS_DIR="./baseEnvs"
FULLSMID="$1"
ENV_FILE="$ENVS_DIR/${FULLSMID}.env"
echo -e "FULLSMID=${FULLSMID}" > $ENV_FILE
echo -e "STAGE_INDIR=/storage1/fs1/cruchagac/Active/$USER/c1in/${FULLSMID}" >> $ENV_FILE
[ ! -d /scratch1/fs1/cruchagac/${USER}/c1in/$FULLSMID ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1in/$FULLSMID
echo -e "INDIR=/input/${FULLSMID}" >> $ENV_FILE
[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/$FULLSMID ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1out/$FULLSMID
echo -e "OUTDIR=/output/${FULLSMID}" >> $ENV_FILE
[ ! -d /storage1/fs1/cruchagac/Active/$USER/c1out/${FULLSMID} ] && mkdir /storage1/fs1/cruchagac/Active/$USER/c1out/${FULLSMID}
echo -e "FINAL_OUTDIR=/final_output/${FULLSMID}" >> $ENV_FILE
echo -e "LOGFILE=/output/${FULLSMID}/${FULLSMID}_runlog.txt" >> $ENV_FILE
echo -e "RUN_TYPE=paddedexome" >> $ENV_FILE
echo -e "RGBASES="$(basename -s .rgfile /storage1/fs1/cruchagac/Active/${USER}/c1in/${FULLSMID}/*.rgfile)"" >> $ENVS_DIR/${FULLSMID}.env
cat ${BASE_ENVS_DIR}/pipelinebase.env >> $ENV_FILE
cat ${BASE_ENVS_DIR}/references.env >> $ENV_FILE
[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/logs/${FULLSMID} ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1out/logs/${FULLSMID}
rsync $ENV_FILE /scratch1/fs1/cruchagac/${USER}/c1out/$FULLSMID/
rsync $ENV_FILE /scratch1/fs1/cruchagac/${USER}/c1out/logs/$FULLSMID/
#if [[ -z $RGBASES ]]; then exit 1; fi
