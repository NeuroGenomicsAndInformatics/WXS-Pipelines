#!bin/bash
SAMPLES_FILE=$1
RUN_NAME=$2
PROJECTS_DIR="/40/AD/AD_Seq_Data/02.-BWA_GATK/02-GRCh38/02-gVCFs"
STAGE_DIR="/storage1/fs1/cruchagac/Active/${USER}/c1in/${RUN_NAME}"
if [ ! -e $STAGE_DIR ]; then mkdir -p ${STAGE_DIR}; fi
CHECK_GVCFS="${STAGE_DIR}/checkgvcf.log"
echo "" > $CHECK_GVCFS
while IFS=$'\t' read -r -a sampArray; do
WXS=${sampArray[0]}
PR=${sampArray[1]}
SM=$(echo ${sampArray[2]} | cut -d. -f1)
if [[ ${WXS^^} == "WES" ]]; then \
rsync -r --include="*/" --include="*vcf*" --exclude="*" $PROJECTS_DIR/${PR}/${SM/MAP_/}* $STAGE_DIR/
elif [[ ${WXS^^} != "WGS" ]]; then echo ${sampArray[@]}  >> $CHECK_GVCFS
fi
done < ${SAMPLES_FILE}
