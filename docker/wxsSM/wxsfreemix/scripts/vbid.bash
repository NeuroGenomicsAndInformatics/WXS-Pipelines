#!/bin/bash
VerifyBamID \
  --BamFile $1 \
  --SVDPrefix /VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat \
  --Reference ${REF_FASTA} \
  --NumThread $LSB_MAX_NUM_PROCESSORS \
  --max-depth 1000 \
