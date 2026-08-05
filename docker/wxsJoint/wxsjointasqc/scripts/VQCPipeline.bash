#!/bin/bash
INT=$1
SNP_RECAL_VCF=$(/scripts/ApplyVQSR.bash ${OUTDIR}/${INT}.joint.vcf.gz SNP ${SNP_RECAL} | tail -n1)
echo -e "\n\n ${SNP_RECAL_VCF} \n\n"
if [[ -f ${INDEL_RECAL} ]]; then
BOTH_RECAL_VCF=$(/scripts/ApplyVQSR.bash ${SNP_RECAL_VCF} INDEL ${INDEL_RECAL} | tail -n1)
else
BOTH_RECAL_VCF=${SNP_RECAL_VCF}
fi
echo -e "\n\n ${BOTH_RECAL_VCF} \n\n"
/scripts/hardfilter.bash ${BOTH_RECAL_VCF}