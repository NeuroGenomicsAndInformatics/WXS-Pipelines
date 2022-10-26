#!/bin/bash
THREADS=$(( LSB_MAX_NUM_PROCESSORS * 2 ))
FQ1=$(head -n1 ${INDIR}/infqfile.txt | cut -d ' ' -f1)
echo $FQ1
RGFILE="${OUTDIR}/${FULLSMID}.$(echo ${FQ1##*/} | cut -d_ -f1 | cut -d. -f1)_$(echo ${FQ1##*/} | cut -d_ -f1 | cut -d. -f2).rgfile"
bwa-mem2 mem -M -t $THREADS -K 10000000 \
  -R $(head -n1 ${RGFILE}) \
  ${REF_FASTA} \
  ${FQ1} \
  ${FQ1/_1.fastq/_2.fastq} \
  | samtools view -b -1 -o ${INDIR}/${FQ1##*/}.aln.bam \
  && rm ${FQ1} && rm ${FQ1/_1.fastq/_2.fastq}
