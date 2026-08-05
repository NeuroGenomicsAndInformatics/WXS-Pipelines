#!/bin/bash
if [[ -n $(find $INDIR -maxdepth 1 -name "*.f*q.gz") ]]; then
  bash /scripts/align_fqs.bash
elif [[ -n $(find $INDIR -maxdepth 1 -name "*.cram") ]]; then
  bash /scripts/align_cram.bash "$1"
elif [[ -n $(find $INDIR -maxdepth 1 -name "*.bam") ]]; then
  bash /scripts/align_bam.bash
else
  exit 74
fi
