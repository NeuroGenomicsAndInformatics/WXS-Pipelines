#!/bin/bash
ENVS_DIR="/scratch1/fs1/cruchagac/matthewj/c1in/envs"
for FULLSMID in $(cat $1) ; do
echo -e "FULLSMID=${FULLSMID}" > $ENVS_DIR/${FULLSMID}.env
echo -e "INDIR=/input/${FULLSMID}" >> $ENVS_DIR/${FULLSMID}.env
mkdir /scratch1/fs1/cruchagac/matthewj/c1out/$FULLSMID
echo -e "OUTDIR=/output/${FULLSMID}" >> $ENVS_DIR/${FULLSMID}.env
echo -e "LOGFILE=/output/${FULLSMID}/log.txt" >> $ENVS_DIR/${FULLSMID}.env
echo -e "RUN_TYPE=paddedexome" >> $ENVS_DIR/${FULLSMID}.env
echo -e "RGBASES="$(basename -s .rgfile /scratch1/fs1/cruchagac/matthewj/c1in/${FULLSMID}/*.rgfile)""
done
