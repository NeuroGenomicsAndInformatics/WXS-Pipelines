#!/bin/bash
JOINT_VCF=$1
SNP_RECAL_TABLE=$2
INDEL_RECAL_TABLE=$3
SNP_RECAL_VCF=$(/scripts/ApplyVQSR_nonpipe.bash ${JOINT_VCF} SNP ${SNP_RECAL_TABLE} ${JOINT_VCF} | tail -n1)
echo -e "\n\n ${SNP_RECAL_VCF} \n\n"
if [[ -s ${INDEL_RECAL_TABLE} ]]; then
BOTH_RECAL_VCF=$(/scripts/ApplyVQSR_nonpipe.bash ${SNP_RECAL_VCF} INDEL ${INDEL_RECAL_TABLE} ${JOINT_VCF} | tail -n1)
else
BOTH_RECAL_VCF=${SNP_RECAL_VCF}
fi
echo -e "\n\n ${BOTH_RECAL_VCF} \n\n"
/scripts/Vfilter1_nonpipe.bash ${BOTH_RECAL_VCF}