#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENVS_DIR="/scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/envs"
[ ! -d $ENVS_DIR ] && mkdir $ENVS_DIR
COHORT="$1"
CHR="$2"
ENV_FILE="$ENVS_DIR/${COHORT}_${CHR}.env"
echo -e "COHORT=${COHORT}" > $ENV_FILE
echo -e "CHR=${CHR}" >> $ENV_FILE
echo -e "INDIR=/storage1/fs1/${STORAGE_USER}/Active/${USER}/c1in/${COHORT}" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR}
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR}/joint_vcfs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR}/joint_vcfs
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR}/dbs ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR}/dbs
echo -e "OUTDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR}" >> $ENV_FILE
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}_${CHR} ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}_${CHR}
echo -e "FINAL_OUTDIR=/storage1/fs1/${STORAGE_USER}/Active/${USER}/c1out/${COHORT}_${CHR}" >> $ENV_FILE
echo -e "JOINT_VCF=${FINAL_OUTDIR}/${COHORT}.${CHR}.wgs.joint.vcf.gz" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}_${CHR} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${COHORT}_${CHR}
cat $SCRIPT_DIR/baseEnvs/references.env >> $ENV_FILE
rsync $ENV_FILE /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${COHORT}_${CHR}
echo $ENV_FILE
