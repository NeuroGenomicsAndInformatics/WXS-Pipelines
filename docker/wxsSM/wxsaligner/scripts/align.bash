#!/bin/bash
if [[ -n $(find $STAGE_INDIR -name "*.cram") ]]; then
  bash /scripts/align_cram.bash "$1"
elif [[ -n $(find $STAGE_INDIR -name "*.bam") ]]; then
  bash /scripts/align_bam.bash
else
  bash /scripts/align_fqs.bash
fi
