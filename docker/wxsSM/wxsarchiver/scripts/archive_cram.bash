#!/bin/bash
samtools view -C -T $REF_FASTA -@ $($LSB_MAX_NUM_PROCESSORS -1) --output-fmt-option archive --output-fmt-option embed_ref -o ${1%.cram}.arc_er.cram $1
