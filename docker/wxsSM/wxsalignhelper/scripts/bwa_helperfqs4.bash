#!/bin/bash
FQ1S=()
for FILE in $(cat ${OUTDIR}/FQ1s.txt); do
FQ1S+=($FILE)
done
FQ1=${FQ1S[$LSB_JOBINDEX-1]}
echo $FQ1
RGFILE="${OUTDIR}/${FULLSMID}.$(echo ${FQ1##*/} | cut -d_ -f1 | cut -d. -f1)_$(echo ${FQ1##*/} | cut -d_ -f1 | cut -d. -f2).rgfile"
if [[ $(wc -c $FQ1 | cut -d' ' -f1) -lt 12000000000 ]]; then
bwa-mem2 mem -M -t $LSB_MAX_NUM_PROCESSORS -K 10000000 \
  -R $(head -n1 ${RGFILE}) \
  ${REF_FASTA} \
  ${FQ1} \
  ${FQ1/_1.fastq/_2.fastq} \
  | samtools sort -@ 4 -m 30G -o ${INDIR}/${FQ1##*/}.bam \
  && rm ${FQ1} && rm ${FQ1/_1.fastq/_2.fastq}
else
bwa-mem2 mem -M -t $LSB_MAX_NUM_PROCESSORS -K 10000000 \
  -R $(head -n1 ${RGFILE}) \
  ${REF_FASTA} \
  ${FQ1} \
  ${FQ1/_1.fastq/_2.fastq} \
  | samtools view -b -1 -o ${INDIR}/${FQ1##*/}.aln.bam \
  && samtools sort -@ 4 -m 30G -o ${INDIR}/${FQ1##*/}.bam ${INDIR}/${FQ1##*/}.aln.bam \
&& rm ${FQ1} && rm ${FQ1/_1.fastq/_2.fastq}
rm ${INDIR}/${FQ1##*/}.aln.bam
fi
