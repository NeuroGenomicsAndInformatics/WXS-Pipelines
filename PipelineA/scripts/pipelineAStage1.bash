#!/bin/bash
source /scripts/pipelineABHelperFunctions.bash
#1. set variables, equivalent to setting the environment in the original pipeline
SAMPLEID=$(echo $FULLSMID | cut -d '^' -f 1)
RGBASE="$(echo ${RGBASES} | cut -d ' ' -f $LSB_JOBINDEX)"
FQ1="$(ls $INDIR/${RGBASE}*1*)"
FQ1EXT="$(echo ${FQ1##*${RGBASE}})"
FQ2="$(echo $INDIR/${RGBASE}${FQ1EXT/1/2})"
echo -e "${FULLSMID}" > ${LOGFILE}
MEM_SPLIT=$((${MEM}/${THREADS}))
echo -e "" > ${OUTDIR}/stage1complete.txt
#TODO Add pipeline B, C, D, and E logic
reportToLog "Starting pipeline A for $RGBASE. Aligning and sorting"
alignSortPairedFQs
#saveToOutputDirectory ${CURRENT_BAM}
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
#saveToOutputDirectory ${CURRENT_BAM}
validateCurrentBam
reportToLog "Validated. Adding to stage1complete list."
#NOTE: This can lead to a race condition on the list. A failed Stage 2 could be caused by the stage1complete.txt file being messed up. Fixing the file manually solves this and Stage 2 can be run normally after.
echo -e "${CURRENT_BAM##*/}" >> ${OUTDIR}/stage1complete.txt
reportToLog "Finished for $RGBASE"
