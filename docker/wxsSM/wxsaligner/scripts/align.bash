#!/bin/bash
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then exit 66
else bash /scripts/gpualign.bash
fi