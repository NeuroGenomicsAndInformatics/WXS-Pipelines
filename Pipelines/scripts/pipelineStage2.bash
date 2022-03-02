#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/pipelineStage2HelperFunctions.bash
stageDataForSample
SAMPLEID=$(echo $FULLSMID | cut -d'^' -f1)
MD_INPUTS=()
sort ${INDIR}/stage1complete.txt | uniq -u > ${OUTDIR}/stage1complete.txt
for BAM in $(cat ${OUTDIR}/stage1complete.txt) ; do
  MD_INPUTS+=("I=${INDIR}/${BAM/CHECK_/}")
done
SAMPLEID_VE=$(echo ${SAMPLEID} | tr "^" "-")
MEM_SPLIT=$((${S2MEM}/${S2THREADS}))
reportToLog "Marking duplicates"
markDuplicates
reportToLog "Marked duplicates. Validating"
validateCurrentBam
reportToLog "Validated. Analyzing depth of coverage."
analyzeDepthOfCoverage
SAMPLE_COVERAGE=$(getAverageDepthOfCoverage)
reportToLog "Average deptho of coverage for ${FULLSMID} is ${SAMPLE_COVERAGE}"
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
if [[ $(wc -c $CURRENT_VCF | cut -d' ' -f1) -lt 40000000 ]]; then
mv ${FINAL_OUTDIR} /final_output/CHECK_${FULLSMID}; cleanUp; exit 8
else rm -R ${STAGE_INDIR}; cleanUp
fi
echo "${SAMPLEID},${SAMPLE_COVERAGE},${SAMPLE_FREEMIX},${SAMPLE_TITV}" > ${FINAL_OUTDIR}/${FULLSMID}_stats.csv
reportToLog "Finished for $FULLSMID."
