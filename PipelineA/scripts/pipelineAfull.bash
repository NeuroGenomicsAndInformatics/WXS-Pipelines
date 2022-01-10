#!/bin/bash
source /scripts/pipelineABHelperFunctions.bash
#1. set variables, equivalent to setting the environment in the original pipeline
SAMPLEID=$(echo $FULLSMID | cut -d '^' -f 1)
RGBASES="$(ls $INDIR/*.rgfile | xargs -I{} basename {} .rgfile)"
MD_INPUTS=()
SAMPLEID_VE=$(echo ${SAMPLEID} | tr "^" "-")
echo -e "${FULLSMID}" > ${LOGFILE}
MEM_SPLIT=$(expr $MEM / $THREADS)
#TODO Add pipeline B, C, D logic
for RGBASE in ${RGBASES[@]}; do
reportToLog "Starting pipeline A for $RGBASE. Aligning and sorting"
alignSortPairedFQs
saveToOutputDirectory ${CURRENT_BAM}
reportToLog "Aligned FASTQs into BAM. Validating"
validateCurrentBam
reportToLog "Validated."
#SAMPLE_FREEMIX=$(getFreeMix)
#reportToLog "FREEMIX for ${FULLSMID} is ${SAMPLE_FREEMIX}"
#if [ ${SAMPLE_FREEMIX} -le 0.03 ]
#then reportToLog "${FULLSMID} is likely contaminated"; exit 3; fi
if [ ${RUN_TYPE} -eq "paddedexome" ]
then
  reportToLog "Intersecting BAM with BED"
  intersectBamWithBed ${CURRENT_BAM} ${REF_PADBED}
  saveToOutputDirectory ${CURRENT_BAM}
  reportToLog "Intersected."
fi
reportToLog "Validating and saving BAM as CRAM"
validateCurrentBam
saveBamAsCram ${CURRENT_BAM}
reportToLog "Saved CRAM. Finished for $RGBASE"
MD_INPUTS+=("I=${CURRENT_BAM}")
done
reportToLog "Marking duplicates"
markDuplicates
reportToLog "Marked duplicates. Validating"
validateCurrentBam
saveToOutputDirectory ${CURRENT_BAM}
reportToLog "Validated. Analyzing depth of coverage."
analyzeDepthOfCoverage
reportToLog "Analyzed. Recalibrating bases"
recalibrateBases
reportToLog "Recalibrated. Validating"
validateCurrentBam
reportToLog "Validated. Calling variants on sample"
callSampleVariants
saveToOutputDirectory ${CURRENT_VCF}
reportToLog "Called. Evaluating variants"
evaluateSampleVariants
reportToLog "Evaluated."
