#!/bin/bash
ENVS_DIR="/scratch1/fs1/cruchagac/matthewj/c1in"
for FULLSMID in $(cat $1) ; do
<<<<<<< HEAD
ENV_FILE="$ENVS_DIR/${FULLSMID}.env"
=======
ENV_FILE="$ENVS_DIR/${FULLSMID}/${FULLSMID}.env"
>>>>>>> cdb7e955411800735c4d4068a5fabd9522b32fb3
echo -e "FULLSMID=${FULLSMID}" > $ENV_FILE
echo -e "INDIR=/input/${FULLSMID}" >> $ENV_FILE
mkdir /scratch1/fs1/cruchagac/matthewj/c1out/$FULLSMID
echo -e "OUTDIR=/output/${FULLSMID}" >> $ENV_FILE
echo -e "LOGFILE=/output/${FULLSMID}/log.txt" >> $ENV_FILE
echo -e "RUN_TYPE=paddedexome" >> $ENV_FILE
echo -e "RGBASES="$(basename -s .rgfile /scratch1/fs1/cruchagac/matthewj/c1in/${FULLSMID}/*.rgfile)"" >> $ENVS_DIR/${FULLSMID}.env
cp $ENV_FILE /scratch1/fs1/cruchagac/matthewj/c1out/$FULLSMID/
done
