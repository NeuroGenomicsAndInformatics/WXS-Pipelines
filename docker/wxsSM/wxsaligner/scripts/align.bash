#!/bin/bash
bash /scripts/stageinfqs.bash
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then
  bash /scripts/cpualign.bash
else
  bash /scripts/gpualign.bash
fi
