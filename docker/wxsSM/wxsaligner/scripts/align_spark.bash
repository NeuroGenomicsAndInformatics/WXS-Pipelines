#!/bin/bash
bash /scripts/stageinfqs.bash
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then
  bash /scripts/cpualign_spark.bash
else
  bash /scripts/gpualign.bash
fi
