#!/bin/bash
for VAR in $(printenv | grep CUDA_VISIBLE_DEVICES); do
export ${VAR/CUDA/NVIDIA}
done
bash /scripts/stageinfqs.bash ${FULLSMID}
pbrun fq2bam \
    --ref ${REF_FASTA} \
    --in-fq-list ${INDIR}/infqfile.txt \
    --out-bam ${INDIR}/${BAM} \
    --out-duplicate-metrics ${METDIR}/${FULLSMID}.dup.metrics.txt \
    --num-gpus 2 \
    --tmp-dir ${TMP_DIR}
