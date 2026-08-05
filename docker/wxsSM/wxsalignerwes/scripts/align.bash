#!/bin/bash
if [[ -n $(find $INDIR -maxdepth 1 -name "*.cram") ]]; then
  export UNWRAP_FASTA="$1"
  bash /scripts/stageincram_cpu.bash $UNWRAP_FASTA
elif [[ -n $(find $INDIR -maxdepth 1 -name "*.bam") ]]; then
  bash /scripts/stageinbam.bash
else
  bash /scripts/stageinfqs.bash
fi
bash /scripts/cpualign.bash