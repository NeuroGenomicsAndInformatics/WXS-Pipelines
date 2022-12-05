#!/bin/bash
if [[ -z $CUDA_VISIBLE_DEVICES ]]; then
  bash /scripts/stageinbam_spark.bash \
  && bash /scripts/cpualign.bash
else
  bash /scripts/stageinbam_gpu.bash \
  && bash /scripts/gpualign.bash
fi
