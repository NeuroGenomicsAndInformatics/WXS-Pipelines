#!/bin/bash
OUTDIR=$2
[ -z $2 ] && OUTDIR=${1%/*}
INFILE=${1##*/}
OUTPUT=${OUTDIR}/${INFILE/%.cram/.arc_er.cram}
samtools view -C -T $REF_FASTA -@ $(( LSB_MAX_NUM_PROCESSORS -1 )) --output-fmt-option archive --output-fmt-option embed_ref -o $OUTPUT $1
