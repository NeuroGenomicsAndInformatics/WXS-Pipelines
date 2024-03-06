#!/bin/bash
if [[ -n $(find $STAGE_INDIR -name "*.cram") ]]; then
  export UNWRAP_FASTA="$1"
  bash /scripts/stageincram_cpu.bash $UNWRAP_FASTA
elif [[ -n $(find $STAGE_INDIR -name "*.bam") ]]; then
  bash /scripts/stageinbam.bash
else
  bash /scripts/stageinfqs.bash
fi
&& bash /scripts/cpualign.bash