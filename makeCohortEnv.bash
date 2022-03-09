#!/bin/bash
ENVS_DIR="/scratch1/fs1/cruchagac/$USER/c1in/envs"
BASE_ENVS_DIR="./baseEnvs"
COHORT="$1"
#SAMPLE_MAP="$2"
INTERVAL="$2"
ENV_FILE="$ENVS_DIR/${COHORT}_${INTERVAL}.env"
echo -e "COHORT=${COHORT}" > $ENV_FILE
echo -e "INTERVAL=${INTERVAL}" >> $ENV_FILE
echo -e "STAGE_INDIR=/storage1/fs1/cruchagac/Active/$USER/c1in/${COHORT}" >> $ENV_FILE
if [ ! -e /scratch1/fs1/cruchagac/$USER/c1in/${COHORT} ]; then mkdir /scratch1/fs1/cruchagac/$USER/c1in/${COHORT}; fi
echo -e "INDIR=/input/${COHORT}" >> $ENV_FILE
if [ ! -e /scratch1/fs1/cruchagac/$USER/c1out/${COHORT}_${INTERVAL} ]; then mkdir /scratch1/fs1/cruchagac/$USER/c1out/${COHORT}_${INTERVAL}; fi
echo -e "OUTDIR=/output/${COHORT}_${INTERVAL}" >> $ENV_FILE
if [ ! -e /storage1/fs1/cruchagac/Active/$USER/c1out/${COHORT}_${INTERVAL} ]; then mkdir /storage1/fs1/cruchagac/Active/$USER/c1out/${COHORT}_${INTERVAL}; fi
echo -e "FINAL_OUTDIR=/final_output/${COHORT}_${INTERVAL}" >> $ENV_FILE
echo -e "LOGFILE=/output/${COHORT}_${INTERVAL}/${COHORT}_${INTERVAL}_runlog.txt" >> $ENV_FILE
#echo -e "SAMPLE_MAP=$2" >> $ENV_FILE
if [[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/logs/${COHORT}_${INTERVAL} ]]; then mkdir /scratch1/fs1/cruchagac/${USER}/c1out/logs/${COHORT}_${INTERVAL}; fi
cat ${BASE_ENVS_DIR}/pipelineBaseJoint.env >> $ENV_FILE
cat ${BASE_ENVS_DIR}/references.env >> $ENV_FILE
rsync $ENV_FILE /scratch1/fs1/cruchagac/$USER/c1out/${COHORT}_${INTERVAL}
