#!/bin/bash
SNP_RECAL_TABLE=$(/scripts/VQSR_SNP_nonpipe.bash | tail -n1)
echo -e "\n\n ${SNP_RECAL_TABLE} \n\n"
INDEL_RECAL_TABLE=$(/scripts/VQSR_INDEL_nonpipe.bash | tail -n1)
echo -e "\n\n ${INDEL_RECAL_TABLE} \n\n"
SNP_RECAL_VCF=$(/scripts/ApplyVQSR_nonpipe.bash ${JOINT_VCF} SNP ${SNP_RECAL_TABLE} | tail -n1)
echo -e "\n\n ${SNP_RECAL_VCF} \n\n"
if [[ -s ${INDEL_RECAL_TABLE} ]]; then
BOTH_RECAL_VCF=$(/scripts/ApplyVQSR_nonpipe.bash ${SNP_RECAL_VCF} INDEL ${INDEL_RECAL_TABLE} | tail -n1)
else
BOTH_RECAL_VCF=${SNP_RECAL_VCF}
fi
echo -e "\n\n ${BOTH_RECAL_VCF} \n\n"