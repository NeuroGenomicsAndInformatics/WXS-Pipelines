#!/bin/bash
source /scripts/pipelineABHelperFunctions.bash
#1. set variables, equivalent to setting the environment in the original pipeline
SAMPLEID=$(echo $FULLSMID | cut -d '^' -f 1)
MD_INPUTS="$(ls $OUTDIR/*.bam | xargs -I{} echo "I={}")"
SAMPLEID_VE=$(echo ${SAMPLEID} | tr "^" "-")
MEM_SPLIT=$((${MEM}/${THREADS}))
#TODO Add pipeline B, C, D logic
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
