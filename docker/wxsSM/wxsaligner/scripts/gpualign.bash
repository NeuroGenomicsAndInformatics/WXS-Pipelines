#!/bin/bash
/usr/local/parabricks/pbrun fq2bam \
    --ref ${REF_FASTA} \
    --in-fq-list ${INDIR}/infqfile.txt \
    --out-bam ${INDIR}/${BAM} \
    --out-duplicate-metrics ${OUTDIR}/${FULLSMID}.dup.metrics.txt \
    --num-gpus 1 \
    --memory-limit 200 \
    --tmp-dir ${TMP_DIR} \
    || (rm -R $INDIR && exit 66)
samtools view -@ $LSB_MAX_NUM_PROCESSORS -C -T ${REF_FASTA} -o $OUTDIR/$CRAM $INDIR/$BAM
rm -R $INDIR
samtools index -@ $LSB_MAX_NUM_PROCESSORS $OUTDIR/$CRAM
