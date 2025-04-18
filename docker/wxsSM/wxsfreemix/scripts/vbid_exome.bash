#!/bin/bash
if [ -z $FULLSMID ]; then FINAL_OUTDIR=${1%/*} && CRAM=${1##*/}; fi
VerifyBamID \
  --BamFile ${FINAL_OUTDIR}/${CRAM} \
  --SVDPrefix /scripts/1000g.phase3.10k.b38.exome.vcf.gz.dat \
  --Reference ${REF_FASTA} \
  --NumThread $LSB_MAX_NUM_PROCESSORS \
  --Output ${FINAL_OUTDIR}/${CRAM}.vbid2 \
  --max-depth 1000
