#!/bin/bash
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then bash /scripts/stageinfqsalign3.bash
else bash gpualign.bash
fi
