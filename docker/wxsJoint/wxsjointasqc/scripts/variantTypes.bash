#!/bin/bash
VCF="$1"

echo "${VCF}" > ${VCF%.*.*}.vcftools.log; zcat ${VCF} | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c >> ${VCF%.*.*}.vcftools.log
