#!/bin/bash
source /scripts/pipelineCommonFunctions.bash
source /scripts/pipelineCHelperFunctions.bash
CRAM="$1"
reportToLog "Starting pipeline C for $CRAM. Getting CRAM"
getCram
reportToLog "Converting CRAM to BAM"
convertCramToBam
reportToLog "CRAM Converted."
#validateCurrentBam
#reportToLog "Validated. Generating rgfiles."
reportToLog "Converting BAM to FASTQs"
revertBamToFastqs
reportToLog "Converted to FASTQs. Making RGfiles and renaming FASTQs"
makeRgfiles
reportToLog "Finished for $CRAM. Adding to Stage 1 workfile."
touch $STAGE1_WORKFILE
echo ${FULLSMID} >> $STAGE1_WORKFILE
