#!/bin/bash
# This script removes monomorphic variants from a vcf file.

INPUT_VCF=$1
bcftools view \
    -e ' INFO/AF = 0 || INFO/AF = 1' \
    -Oz -o ${INPUT_VCF%.*.*}.nomono.vcf.gz \
    $INPUT_VCF

#bcftools view -e 'INFO/AF = 0 || INFO/AF = 1' -Oz -o test.chr9:93910711-98853379.vcf.gz gcad.preview.r5.wgs.58507.GLnexus.2024.11.03.genotypes.chr9:93910711-98853379.ALL.vcf.bgz

#bcftools view -e ' COUNT(GT="RR") = INFO/AN / 2 || COUNT( GT="AA") = INFO/AN / 2' -Oz -o ${INPUT_VCF%.*.*}.nomono.vcf.gz $INPUT_VCF