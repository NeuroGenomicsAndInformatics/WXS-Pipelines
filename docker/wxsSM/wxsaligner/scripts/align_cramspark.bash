#!/bin/bash
export UNWRAP_FASTA="$1"
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then
  bash /scripts/stageincram_spark.bash $UNWRAP_FASTA \
  && bash /scripts/cpualign_spark.bash
else
  bash /scripts/stageincram_gpu.bash $UNWRAP_FASTA \
  && bash /scripts/gpualign.bash
fi