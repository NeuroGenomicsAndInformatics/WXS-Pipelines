#!/bin/bash
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then bash /scripts/stageinfqsalign3.bash
else bash /scripts/gpualign.bash
fi
