#!/bin/bash
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then
  bash stageinfqsalign3spark.bash
else
  bash /scripts/gpualign.bash
fi
