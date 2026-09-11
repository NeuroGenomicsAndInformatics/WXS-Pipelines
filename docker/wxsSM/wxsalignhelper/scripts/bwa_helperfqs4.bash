#!/bin/bash
THREADS=$(( LSB_MAX_NUM_PROCESSORS * 2 ))
FQ1=$(head -n $LSB_JOBINDEX ${INDIR}/infqfile.txt | tail -n1 | cut -d ' ' -f1)
echo $FQ1
#RGFILE="${OUTDIR}/${FULLSMID}.$(echo ${FQ1##*/} | cut -d_ -f1 | cut -d. -f1)_$(echo ${FQ1##*/} | cut -d_ -f1 | cut -d. -f2).rgfile"
if [[ $(wc -c $FQ1 | cut -d' ' -f1) -lt 5000000000 ]]; then
bwa-mem2 mem -M -t $THREADS -K 10000000 \
  -R $(head -n $LSB_JOBINDEX ${INDIR}/infqfile.txt | tail -n1 | cut -d ' ' -f3) \
  ${REF_FASTA} \
  ${FQ1} \
  ${FQ1/1.f/2.f} \
  | ${GATK4261mod} \
  --java-options "-Xmx30g -XX:ParallelGCThreads=2" \
  SortSam  \
  -I /dev/stdin \
  -O ${FQ1}.bam \
  -R ${REF_FASTA} \
  -SO coordinate \
  --MAX_RECORDS_IN_RAM 1000000 \
  --CREATE_INDEX true \
  --TMP_DIR $TMP_DIR
else
bwa-mem2 mem -M -t $THREADS -K 10000000 \
  -R $(head -n $LSB_JOBINDEX ${INDIR}/infqfile.txt | tail -n1 | cut -d ' ' -f3) \
  ${REF_FASTA} \
  ${FQ1} \
  ${FQ1/1.f/2.f} \
  | samtools view -b -o ${FQ1}.aln.bam \
  && ${GATK4261mod} \
  --java-options "-Xmx30g -XX:ParallelGCThreads=2" \
  SortSam  \
  -I ${FQ1}.aln.bam \
  -O ${FQ1}.bam \
  -R ${REF_FASTA} \
  -SO coordinate \
  --MAX_RECORDS_IN_RAM 1000000 \
  --CREATE_INDEX true \
  --TMP_DIR $TMP_DIR \
&& rm ${FQ1}.aln.bam
fi
[[ -s ${FQ1}.bai ]] && rm ${INDIR}/${LSB_JOBINDEX}.lock
