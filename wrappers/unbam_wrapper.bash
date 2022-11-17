#!/bin/bash
BAM=$1

[ -f ${BAM}.bai ] || samtools index $BAM

gatk --java-options "-Xmx40g -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  RevertSam \
    -I ${BAM} \
    -O /dev/stdout \
    -SO queryname \
    --COMPRESSION_LEVEL 0 \
    --VALIDATION_STRINGENCY SILENT
| gatk --java-options "-Xmx40g -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SamToFastq \
    -I /dev/stdin \
    --COMPRESS_OUTPUTS_PER_RG true \
    --OUTPUT_PER_RG true \
    --OUTPUT_DIR ${BAM%/*} \
    -RG_TAG PU \
    --VALIDATION_STRINGENCY SILENT
