#!/bin/bash
bash /scripts/stageinfqs.bash ${FULLSMID}
pbrun fq2bam \
    --ref ${REF_FASTA} \
    --in-fq-list ${INDIR}/infqfile.txt \
    --out-bam ${INDIR}/${BAM} \
    --out-duplicate-metrics ${METDIR}/${FULLSMID}.dup.metrics.txt \
    --num-gpus 1 \
    --tmp-dir ${TMP_DIR}
samtools view -@ $LSB_MAX_NUM_PROCESSORS -C -T ${REF_FASTA} -o $OUTDIR/$CRAM $INDIR/$BAM
samtools index -@ $LSB_MAX_NUM_PROCESSORS $OUTDIR/$CRAM && rm $INDIR
