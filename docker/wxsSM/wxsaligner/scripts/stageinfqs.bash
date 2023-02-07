#!/bin/bash
rsync -rL $STAGE_INDIR/ $INDIR
sleep 10
INFQ_FILE=${INDIR}/infqfile.txt
echo -n "" > $INFQ_FILE
for FQ in $(find $INDIR -name "*_1.f*q.gz"); do
SM=$(echo $FULLSMID | cut -d^ -f1)
BARCODE=$(echo $FULLSMID | cut -d^ -f2)
PROJECT=$(echo $FULLSMID | cut -d^ -f3)
FLOWCELL=$(zcat ${FQ} | head -n1 | cut -d: -f3)
LANE=$(zcat ${FQ} | head -n1 | cut -d: -f4)
FLOWLANE="${FLOWCELL}.${LANE}"
[[ -z $FLOWCELL ]] && FLOWLANE=$(echo ${FQ##*/} | rev | cut -d_ -f2- | rev)
echo "@RG\tID:${FLOWLANE}\tPL:illumina\tPU:${FLOWLANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" > ${OUTDIR}/${FULLSMID}.${FLOWLANE}.rgfile
echo "${FQ} ${FQ/_1.f/_2.f} @RG\tID:${FLOWLANE}\tPL:illumina\tPU:${FLOWLANE}.${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" >> ${INFQ_FILE}
done