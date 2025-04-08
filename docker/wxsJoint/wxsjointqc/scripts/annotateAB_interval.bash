#!/bin/bash
# This script adds standard annotations to a vcf file.
INPUT_DIR=$1
INTERVAL=$2
[[ -z $INTERVAL ]] && INTERVAL=${LSB_JOBINDEX}
INPUT_VCF=$(find ${INPUT_DIR} -name "*.${INTERVAL}.setmiss.vcf.gz")
INT_LISTS=($(ls ${INPUT_DIR}/intlists | grep scattered.interval_list))
INT_LIST="${INPUT_DIR}/intlists/${INT_LISTS[${INTERVAL}-1]}"
OUTPUT_VCF="${INPUT_VCF%.*.*}.annAB.vcf.gz"

java -Xmx18g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=/tmp \
	-jar /ref/GATK360.jar \
	-T VariantAnnotator \
	-A AlleleBalance \
	-R ${REF_FASTA} \
	-V ${INPUT_VCF} \
	-o ${OUTPUT_VCF} 

SUCCESS=$?

[[ $SUCCESS -eq 0 ]] && rm ${INPUT_VCF}*

exit $SUCCESS