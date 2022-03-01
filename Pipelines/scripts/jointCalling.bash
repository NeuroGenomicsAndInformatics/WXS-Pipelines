#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/jointCallingHelperFunctions.bash
#1. set variables, equivalent to setting the environment in the original pipeline
MEM_SPLIT=$((${S3MEM}/${S3THREADS}))
reportToLog "Staging Data"
stageDataForCOHORT
reportToLog "Staged. Building Sample Map"
makeSampleMap
reportToLog "Building GenomicDB"
buildGenomicDB
reportToLog "Built. Joint Calling"
jointCallCohort
#saveToOutputDirectory ${CURRENT_BAM}
reportToLog "Called. Transferring files to storage"
transferOutputFilesToStorage
cleanUp
reportToLog "Finished for $COHORT."
