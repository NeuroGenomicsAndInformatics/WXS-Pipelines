#!/bin/bash
export NVIDIA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES
for VAR in $(cat /scratch1/fs1/cruchagac/$USER/c1in/envs/${1}.env); do
export $VAR;
done
bash /scripts/stageinfqs.bash ${FULLSMID}
pbrun fq2bam \
    --ref ${REF_FASTA} \
    --in-fq-list ${INDIR}/infqfile.txt \
    --out-bam ${INDIR}/${BAM} \
    --out-duplicate-metrics ${METDIR}/${FULLSMID}.dup.metrics.txt \
    --num-gpus 1 \
    --tmp-dir ${TMP_DIR} \
    || exit 66
samtools view -@ $LSB_MAX_NUM_PROCESSORS -C -T ${REF_FASTA} -o $OUTDIR/$CRAM $INDIR/$BAM
samtools index -@ $LSB_MAX_NUM_PROCESSORS $OUTDIR/$CRAM && rm $INDIR
