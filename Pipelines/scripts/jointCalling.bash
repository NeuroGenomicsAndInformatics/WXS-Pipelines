#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/jointCallingHelperFunctions.bash
MEM_SPLIT=$((${S3MEM}/${S3THREADS}))
reportToLog "Staging Data"
stageDataForCOHORT
reportToLog "Staged. Building Sample Map"
makeSampleMap
reportToLog "Building GenomicDB"
buildGenomicDB
reportToLog "Built. Joint Calling"
jointCallCohort
reportToLog "Called. Transferring files to storage"
transferOutputFilesToStorage
#cleanUp
reportToLog "Finished for $COHORT."
