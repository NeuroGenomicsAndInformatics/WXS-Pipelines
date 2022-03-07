#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/jointCallingHelperFunctions.bash
#1. set variables, equivalent to setting the environment in the original pipeline
MEM_SPLIT=$((${S3MEM}/${S3THREADS}))
reportToLog "Staging Data"
stageDataForCOHORT
reportToLog "Staged. Building Sample Map"
#makeSampleMap
echo -n "" > ${INDIR}/sampmap_clean.txt
cat ${INDIR}/sampmap.txt | while read LINE; do
  SAMPLE=($LINE)
  echo -n ${SAMPLE[0]//^/.} >> ${INDIR}/sampmap_clean.txt
  echo -e -n "\t" >> ${INDIR}/sampmap_clean.txt
  echo ${SAMPLE[1]//^/.} >> ${INDIR}/sampmap_clean.txt
done
SAMPLE_MAP=${INDIR}/sampmap_clean.txt
reportToLog "Building GenomicDB"
buildGenomicDB
reportToLog "Built. Joint Calling"
jointCallCohort
#saveToOutputDirectory ${CURRENT_BAM}
reportToLog "Called. Transferring files to storage"
transferOutputFilesToStorage
#cleanUp
reportToLog "Finished for $COHORT."
