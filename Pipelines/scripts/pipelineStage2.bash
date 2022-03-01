#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/pipelineStage2HelperFunctions.bash
stageDataForSample
SAMPLEID=$(echo $FULLSMID | cut -d'^' -f1)
MD_INPUTS=()
sort ${INDIR}/stage1complete.txt | uniq -u > ${OUTDIR}/stage1complete.txt
for BAM in $(cat ${OUTDIR}/stage1complete.txt) ; do
  MD_INPUTS+=("I=${OUTDIR}/${BAM/CHECK_/}")
done
SAMPLEID_VE=$(echo ${SAMPLEID} | tr "^" "-")
MEM_SPLIT=$((${S2MEM}/${S2THREADS}))
reportToLog "Marking duplicates"
markDuplicates
reportToLog "Marked duplicates. Validating"
validateCurrentBam
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
reportToLog "Called. Evaluating variants"
evaluateSampleVariants
reportToLog "Evaluated."
SAMPLE_TITV=$(getTitvRatio)
reportToLog "TITV for ${FULLSMID} is ${SAMPLE_TITV}"
reportToLog "Transferring output files to ${FINAL_OUTDIR}"
transferOutputFilesToStorage
cleanUp
if [[ $(wc -c $CURRENT_VCF | cut -d' ' -f1) -lt 40000000 ]]; then
mv ${FINAL_OUTDIR} /final_output/CHECK_${FULLSMID}; exit 8
else rm -R ${STAGE_INDIR}
fi
reportToLog "Finished for $FULLSMID."
