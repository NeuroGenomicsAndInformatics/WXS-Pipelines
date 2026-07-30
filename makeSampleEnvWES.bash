#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENVS_DIR="/storage1/fs1/${STORAGE_USER}/Active/${USER}/c1in/envs"
[ ! -d $ENVS_DIR ] && mkdir $ENVS_DIR
BASE_ENVS_DIR="${SCRIPT_DIR}/baseEnvs"
FULLSMID="$1"
ENV_FILE="$ENVS_DIR/${FULLSMID}.env"
echo -e "ENV_FILE=${ENV_FILE}" > $ENV_FILE
echo -e "FULLSMID=${FULLSMID}" >> $ENV_FILE
echo -e "INDIR=/storage1/fs1/${STORAGE_USER}/Active/$USER/c1in/${FULLSMID}" >> $ENV_FILE
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/$USER/c1in/${FULLSMID} ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/$USER/c1in/${FULLSMID}
echo -e "OUTDIR=/storage1/fs1/${STORAGE_USER}/Active/$USER/c1out/${FULLSMID}" >> $ENV_FILE
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/$USER/c1out/${FULLSMID} ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/$USER/c1out/${FULLSMID}
echo -e "LOGDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}" >> $ENV_FILE
echo -e "REF_DIR=${REF_DIR}" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/tmp/${FULLSMID} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/tmp/${FULLSMID}
echo -e "TMP_DIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/tmp/${FULLSMID}" >> $ENV_FILE
echo -e "RUN_TYPE=padded_exome" >> $ENV_FILE
echo -e "BAM=${FULLSMID}.aln.srt.mrk.bam" >> $ENV_FILE
echo -e "CRAM=${FULLSMID}.aln.srt.mrk.cram" >> $ENV_FILE
echo -e "GVCF=${FULLSMID}.snp.indel.g.vcf.gz" >> $ENV_FILE
cat ${BASE_ENVS_DIR}/references_2_1.env >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}
