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
[ ! -d /scratch1/fs1/cruchagac/$USER/c1in/${COHORT} ] && mkdir /scratch1/fs1/cruchagac/$USER/c1in/${COHORT}
echo -e "INDIR=/input/${COHORT}" >> $ENV_FILE
[ ! -d /scratch1/fs1/cruchagac/$USER/c1out/${COHORT}_${INTERVAL} ] && mkdir /scratch1/fs1/cruchagac/$USER/c1out/${COHORT}_${INTERVAL}
echo -e "OUTDIR=/output/${COHORT}_${INTERVAL}" >> $ENV_FILE
[ ! -d /storage1/fs1/cruchagac/Active/$USER/c1out/${COHORT}_${INTERVAL} ] && mkdir /storage1/fs1/cruchagac/Active/$USER/c1out/${COHORT}_${INTERVAL}
echo -e "FINAL_OUTDIR=/final_output/${COHORT}_${INTERVAL}" >> $ENV_FILE
echo -e "LOGFILE=/output/${COHORT}_${INTERVAL}/${COHORT}_${INTERVAL}_runlog.txt" >> $ENV_FILE
#echo -e "SAMPLE_MAP=$2" >> $ENV_FILE
[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/logs/${COHORT}_${INTERVAL} ] && mkdir /scratch1/fs1/cruchagac/${USER}/c1out/logs/${COHORT}_${INTERVAL}
cat ${BASE_ENVS_DIR}/pipelineBaseJoint.env >> $ENV_FILE
cat ${BASE_ENVS_DIR}/references.env >> $ENV_FILE
rsync $ENV_FILE /scratch1/fs1/cruchagac/$USER/c1out/${COHORT}_${INTERVAL}
