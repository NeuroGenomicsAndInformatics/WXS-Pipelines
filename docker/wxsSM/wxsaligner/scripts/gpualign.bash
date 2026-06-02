#!/bin/bash
rsync -rL $INDIR/ $TMP_DIR/
/usr/local/parabricks/pbrun fq2bam \
    --ref ${REF_FASTA} \
    --in-fq-list ${TMP_DIR}/infqfile.txt \
    --out-bam ${TMP_DIR}/${BAM} \
    --out-duplicate-metrics ${OUTDIR}/${FULLSMID}.dup.metrics.txt \
    --num-gpus 1 \
    --memory-limit 200 \
    --tmp-dir ${TMP_DIR} \
    || (rm -R ${TMP_DIR} && exit 66)
samtools view -@ $LSB_MAX_NUM_PROCESSORS -C -T ${REF_FASTA} -o $OUTDIR/$CRAM ${TMP_DIR}/${FULLSMID}/$BAM
rm -R ${TMP_DIR}
samtools index -@ $LSB_MAX_NUM_PROCESSORS $OUTDIR/$CRAM
