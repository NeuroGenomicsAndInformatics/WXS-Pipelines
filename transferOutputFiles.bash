#!/bin/bash
OUTDIR=/scratch1/fs1/cruchagac/$USER/c1out
FINAL_OUTDIR=/storage1/fs1/cruchagac/Active/$USER/c1out
for FULLSMID in $(cat $1); do
mkdir ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*.env ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*.cram ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*_GATK* ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*vcf* ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*verifybam.selfSM ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*_exome_* ${FINAL_OUTDIR}/${FULLSMID}
done
