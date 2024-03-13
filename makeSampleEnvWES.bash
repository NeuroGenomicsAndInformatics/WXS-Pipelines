#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
ENVS_DIR="/scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/envs"
[ ! -d $ENVS_DIR ] && mkdir $ENVS_DIR
BASE_ENVS_DIR="${SCRIPT_DIR}/baseEnvs"
FULLSMID="$1"
ENV_FILE="$ENVS_DIR/${FULLSMID}.env"
echo -e "ENV_FILE=${ENV_FILE}" > $ENV_FILE
echo -e "FULLSMID=${FULLSMID}" >> $ENV_FILE
echo -e "STAGE_INDIR=/storage1/fs1/${STORAGE_USER}/Active/$USER/c1in/${FULLSMID}" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/$FULLSMID ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/$FULLSMID
echo -e "INDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1in/${FULLSMID}" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/$FULLSMID ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/$FULLSMID
echo -e "OUTDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${FULLSMID}" >> $ENV_FILE
[ ! -d /storage1/fs1/${STORAGE_USER}/Active/$USER/c1out/${FULLSMID} ] && mkdir /storage1/fs1/${STORAGE_USER}/Active/$USER/c1out/${FULLSMID}
echo -e "FINAL_OUTDIR=/storage1/fs1/${STORAGE_USER}/Active/$USER/c1out/${FULLSMID}" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/$FULLSMID ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/$FULLSMID
echo -e "LOGDIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}" >> $ENV_FILE
echo -e "REF_DIR=${REF_DIR}" >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/tmp ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/tmp
echo -e "TMP_DIR=/scratch1/fs1/${SCRATCH_USER}/${USER}/tmp" >> $ENV_FILE
echo -e "RUN_TYPE=padded_exome" >> $ENV_FILE
echo -e "BAM=${FULLSMID}.aln.srt.mrk.bam" >> $ENV_FILE
echo -e "CRAM=${FULLSMID}.aln.srt.mrk.cram" >> $ENV_FILE
echo -e "GVCF=${FULLSMID}.snp.indel.g.vcf.gz" >> $ENV_FILE
cat ${BASE_ENVS_DIR}/references_2_0.env >> $ENV_FILE
[ ! -d /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID} ] && mkdir /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}
rsync $ENV_FILE /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/${FULLSMID}/
rsync $ENV_FILE /scratch1/fs1/${SCRATCH_USER}/${USER}/c1out/logs/${FULLSMID}/
