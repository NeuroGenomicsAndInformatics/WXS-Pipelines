#!bin/bash
SAMP_MAP="$1"
WORKFILE="storage1/fs1/cruchagac/Active/${USER}/workfile$3.txt"
if [[ ! -f $WORKFILE ]]; then echo -n "" > $WORKFILE; fi
FILES_DIR=${SAMP_MAP%/*}
{
  IFS='\n' read
  while IFS='\n' read -r LINE; do
READNUM=$(echo $FILE | cut -d "," -f 5)
if [ "${READNUM}" = "1" ]; then
FLOWCELL=$(echo $FILE | cut -d "," -f 2)
LANE=$(echo $FILE | cut -d "," -f 4)
SM=$(echo $(echo $FILE | cut -d "," -f 6) | cut -d "-" -f 2)
BARCODE=$(echo $(echo $FILE | cut -d "," -f 6) | cut -d "-" -f 3)
PROJECT="$2"
PROJECTNAME="$3"
FULLSM=$(echo ${SM}.${PROJECTNAME}\^${BARCODE}\^${PROJECT})
RGBASE="${FULLSM}.${FLOWCELL}^${LANE}"
STAGE_DIR="/storage1/fs1/cruchagac/Active/${USER}/c1in/${FULLSM}"
if [ ! -e $STAGE_DIR ]; then mkdir ${STAGE_DIR}; echo $FULLSM >> $WORKFILE; fi
RGFILE="$STAGE_DIR/${RGBASE}.rgfile"
touch $RGFILE
echo "@RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}.${PROJECTNAME}\tDS:${FULLSM}" > $RGFILE
#rsync ${RGFILE} ${STAGE_DIR}/${RGBASE}.rgfile
FQ1FILE="$(echo $FILE | cut -d "," -f 1)"
ln -s ${FILES_DIR}/${FQ1FILE} ${STAGE_DIR}/${RGBASE}.r1.fq.gz
ln -s ${FILES_DIR}/${FQ1FILE%_*}_R2.fastq.gz ${STAGE_DIR}/${RGBASE}.r2.fq.gz
fi
done
} < $SAMP_MAP
