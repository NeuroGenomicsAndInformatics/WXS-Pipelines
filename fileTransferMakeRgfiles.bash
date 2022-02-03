#!bin/bash
SAMP_MAP="$1"
FILES_DIR=${SAMP_MAP%/*}
for FILE in $(cat "$1"); do
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
STAGE_DIR="/scratch1/fs1/cruchagac/$USER/c1in/${FULLSM}"
if [ ! -e $STAGE_DIR ]; then mkdir ${STAGE_DIR}; fi
RGFILE="$STAGE_DIR/${RGBASE}.rgfile"
touch $RGFILE
echo "@RG\tID:${FLOWCELL}:${LANE}\tPL:illumina\tPU:${FLOWCELL}:${LANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}.${PROJECTNAME}\tDS:${FULLSM}" > $RGFILE
#rsync ${RGFILE} ${STAGE_DIR}/${RGBASE}.rgfile
FQ1FILE="$(echo $FILE | cut -d "," -f 1)"
rsync ${FILES_DIR}/${FQ1FILE} ${STAGE_DIR}/${RGBASE}.r1.fq.gz
rsync ${FILES_DIR}/${FQ1FILE%_*}_R2.fastq.gz ${STAGE_DIR}/${RGBASE}.r2.fq.gz
fi
done
