#!/bin/bash
source /scripts/pipelineABHelperFunctions.bash
#1. set variables, equivalent to setting the environment in the original pipeline
SAMPLEID=$(echo $FULLSMID | cut -d '^' -f 1)
RGBASE="$(echo ${RGBASES} | cut -d ' ' -f $LSB_JOBINDEX)"
echo -e "${FULLSMID}" > ${LOGFILE}
MEM_SPLIT=$((${MEM}/${THREADS}))
#TODO Add pipeline B, C, D logic
reportToLog "Starting pipeline A for $RGBASE. Aligning and sorting"
alignSortPairedFQs
reportToLog "Aligned FASTQs into BAM. Validating"
validateCurrentBam
reportToLog "Validated."
#SAMPLE_FREEMIX=$(getFreeMix)
#reportToLog "FREEMIX for ${FULLSMID} is ${SAMPLE_FREEMIX}"
#if [ ${SAMPLE_FREEMIX} -le 0.03 ]
#then reportToLog "${FULLSMID} is likely contaminated"; exit 3; fi
if [ "${RUN_TYPE}" = "paddedexome" ]
then
  reportToLog "Intersecting BAM with BED"
  intersectBamWithBed ${CURRENT_BAM} ${REF_PADBED}
  reportToLog "Intersected."
fi
reportToLog "Validating and saving BAM as CRAM"
validateCurrentBam
saveToOutputDirectory ${CURRENT_BAM}
saveBamAsCram ${CURRENT_BAM}
reportToLog "Saved CRAM. Finished for $RGBASE"
