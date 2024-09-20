#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENVS_DIR="/scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/envs"
[ ! -d $ENVS_DIR ] && mkdir $ENVS_DIR
COHORT="$1"
ENV_FILE="$ENVS_DIR/${COHORT}.env"
echo -e "COHORT=${COHORT}" > $ENV_FILE
echo -e "INDIR=/storage1/fs1/${STORAGE_USER}/Active/${USER}/c1in/${COHORT}" >> $ENV_FILE
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}_work ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}_work
echo -e "OUTDIR=/storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}_work" >> $ENV_FILE
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT} ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}
echo -e "FINAL_OUTDIR=/storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}
cat $SCRIPT_DIR/baseEnvs/references_2_0.env >> $ENV_FILE
rsync $ENV_FILE /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}
echo $ENV_FILE
