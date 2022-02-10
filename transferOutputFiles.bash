#!/bin/bash
OUTDIR=/scratch1/fs1/cruchagac/matthewj/c1out
FINAL_OUTDIR=/storage1/fs1/cruchagac/Active/matthewj/c1out
FULLSMID="$1"
mkdir ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*.env ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*.cram ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*_GATK* ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*vcf* ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*verifybam.selfSM ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*_exome_* ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*.pdf ${FINAL_OUTDIR}/${FULLSMID}
rsync ${OUTDIR}/${FULLSMID}/*.recal.table* ${FINAL_OUTDIR}/${FULLSMID}
#rm -R ${OUTDIR}/${FULLSMID}
