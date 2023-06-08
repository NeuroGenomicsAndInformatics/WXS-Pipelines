#!/bin/bash
CRAM=$1
SEX_FILE="${CRAM}.sex.csv"

## Get sex information
echo "CRAM,X coverage,Y coverage,ratio,sex call" > $SEX_FILE
XCOV=$(samtools coverage --reference $REF_FASTA -r chrX ${CRAM} | tail -n1 | cut -f6)
YCOV=$(samtools coverage --reference $REF_FASTA -r chrY ${CRAM} | tail -n1 | cut -f6)
SEX_RATIO=$(bc -l <<<"$XCOV/$YCOV")
SEX_CALL="NA"
if (( $(echo "$SEX_RATIO > 10" | bc -l) )); then SEX_CALL=F; fi
if (( $(echo "$SEX_RATIO < 10" | bc -l) )); then SEX_CALL=M; fi
echo "${CRAM},${XCOV},${YCOV},${SEX_RATIO}" >> $SEX_FILE