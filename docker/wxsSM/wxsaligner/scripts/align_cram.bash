#!/bin/bash
export UNWRAP_FASTA="$1"
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then
  bash /scripts/stageincram_cpu.bash $UNWRAP_FASTA \
  && bash /scripts/cpualign.bash
else
  bash /scripts/stageincram_cpu.bash $UNWRAP_FASTA \
  && bash /scripts/gpualign.bash
fi
