#!/bin/bash
CRAM=$1
UNWRAP_FASTA=$2

[ -f ${CRAM}.crai ] || samtools index $CRAM

gatk --java-options "-Xmx40g -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  RevertSam \
    -I ${CRAM} \
    -O /dev/stdout \
    -R ${UNWRAP_FASTA} \
    -SO queryname \
    --COMPRESSION_LEVEL 0 \
    --VALIDATION_STRINGENCY SILENT
| gatk --java-options "-Xmx40g -XX:ParallelGCThreads=1 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
  SamToFastq \
    -I /dev/stdin \
    --COMPRESS_OUTPUTS_PER_RG true \
    --OUTPUT_PER_RG true \
    --OUTPUT_DIR ${CRAM%/*} \
    -RG_TAG PU \
    --VALIDATION_STRINGENCY SILENT
