#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/pipelineStage2HelperFunctions.bash
#1. set variables, equivalent to setting the environment in the original pipeline
SAMPLEID=$(echo $FULLSMID | cut -d '^' -f 1)
MD_INPUTS=()
for BAM in $(cat ${OUTDIR}/stage1complete.txt) ; do \
MD_INPUTS+=("I=${OUTDIR}/${BAM}")
done
SAMPLEID_VE=$(echo ${SAMPLEID} | tr "^" "-")
MEM_SPLIT=$((${MEM}/${THREADS}))
#TODO Add pipeline B, C, D logic
reportToLog "Marking duplicates"
markDuplicates
reportToLog "Marked duplicates. Validating"
validateCurrentBam
#saveToOutputDirectory ${CURRENT_BAM}
reportToLog "Validated. Analyzing depth of coverage."
analyzeDepthOfCoverage
reportToLog "Analyzed. Recalibrating bases"
recalibrateBases
reportToLog "Recalibrated. Validating"
validateCurrentBam
reportToLog "Validated. Getting FREEMIX."
SAMPLE_FREEMIX=$(getFreeMix)
reportToLog "FREEMIX for ${FULLSMID} is ${SAMPLE_FREEMIX}"
#if [ ${SAMPLE_FREEMIX} -le 0.03 ]
#then reportToLog "${FULLSMID} is likely contaminated"; exit 3; fi
reportToLog "Calling variants on sample"
callSampleVariants
#saveToOutputDirectory ${CURRENT_VCF}
reportToLog "Called. Evaluating variants"
evaluateSampleVariants
reportToLog "Evaluated."
SAMPLE_TITV=$(getTitvRatio)
reportToLog "TITV for ${FULLSMID} is ${SAMPLE_TITV}"
reportToLog "Transferring output files to $FINAL_OUTDIR"
transferOutputFilesToStorage
reportToLog "Finished for $FULLSMID."
