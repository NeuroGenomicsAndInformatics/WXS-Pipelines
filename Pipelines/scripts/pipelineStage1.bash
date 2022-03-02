#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/pipelineStage1HelperFunctions.bash
SAMPLEID=$(echo $FULLSMID | cut -d'^' -f1)
RG_ARRAY=($RGBASES)
RGBASE="${RG_ARRAY[$LSB_JOBINDEX-1]}"
MEM_SPLIT=$((${S1MEM}/${S1THREADS}))
reportToLog "Starting pipeline for $RGBASE. Staging Data to scratch1."
stageDataForRGBASE
FQ1="$(find $INDIR -name "${RGBASE}*1*" -print)"
FQ1EXT="$(echo ${FQ1##*${RGBASE}})"
if [[ -e $FQ1 ]]; then FQ2="$(echo $INDIR/${RGBASE}${FQ1EXT/1/2})"; else unset FQ2; fi
FQI="$(find $INDIR -name "${RGBASE}*f*q*" -print)"
touch ${INDIR}/stage1complete.txt
reportToLog "Aligning and sorting"
if [[ -e $FQ1 ]] && [[ -e $FQ2 ]]; then
if [[ $(wc -c $FQ1 | cut -d' ' -f1) -lt 14000000000 ]]; then alignSortPairedFQs || alignSortPairedHugeFQs; else alignSortPairedHugeFQs; fi
elif [[ ! -e $FQ1 ]] && [[ ! -e $FQ2 ]] && [[ -e $FQI ]]; then
alignSortInterleavedFQs || alignSortHugeInterleavedFQs
else
reportToLog "Check input files. $FQ1 $FQ2 $FQI"; exit 3;
fi
reportToLog "Aligned FASTQs into BAM. Validating"
validateCurrentBam
reportToLog "Validating and saving BAM as CRAM"
saveBamAsCram ${CURRENT_BAM}
if [ "${RUN_TYPE}" = "paddedexome" ]
then
  reportToLog "Intersecting BAM with BED"
  intersectBamWithBed ${CURRENT_BAM} ${REF_PADBED}
  reportToLog "Intersected."
fi
validateCurrentBam
reportToLog "Validated. Adding to stage1complete list."
#NOTE: This can lead to a race condition on the list. A failed Stage 2 could be caused by the stage1complete.txt file being messed up. Fixing the file manually solves this and Stage 2 can be run normally after.
if [[ $(wc -c $CURRENT_BAM | cut -d' ' -f1) -gt 80000000 ]]; then
echo -e "${CURRENT_BAM##*/}" >> ${INDIR}/stage1complete.txt
rm ${STAGE_INDIR}/${RGBASE}*
RETURN_VAL=0
else
echo -e "CHECK_${CURRENT_BAM##*/}" >> ${INDIR}/stage1complete.txt
RETURN_VAL=5
touch ${STAGE_INDIR}/stage2.lock
fi
saveToStageDirectory ${CURRENT_BAM}
saveToStageDirectory ${INDIR}/stage1complete.txt
saveToStageDirectory ${OUTDIR}/${RGBASE}.cram
cleanUp
reportToLog "Finished for $RGBASE"
return RETURN_VAL
