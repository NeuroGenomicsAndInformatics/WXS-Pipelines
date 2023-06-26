#!/bin/bash
CRAM=$1
SEX_FILE="${CRAM}.sex.csv"

## Get sex information
XCOV=$(samtools coverage --reference $REF_FASTA -r chrX ${CRAM} | tail -n1 | cut -f6)
YCOV=$(samtools coverage --reference $REF_FASTA -r chrY ${CRAM} | tail -n1 | cut -f6)
SEX_RATIO=$(bc -l <<<"$XCOV/$YCOV")
SEX_CALL="0"
if (( $(echo "$SEX_RATIO > 15" | bc -l) )); then SEX_CALL=2; fi
if (( $(echo "$SEX_RATIO < 5" | bc -l) )); then SEX_CALL=1; fi
echo "CRAM,X coverage,Y coverage,sex ratio,sex call" > $SEX_FILE
echo "${CRAM##*/},${XCOV},${YCOV},${SEX_RATIO},${SEX_CALL}" >> $SEX_FILE