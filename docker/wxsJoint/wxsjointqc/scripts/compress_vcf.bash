#!/bin/bash

INPUT_DIR=$1
INPUT_VCF=$(find $INPUT_DIR -name "*.vcf.gz")

bcftools view --threads 4 -W=tbi -Oz7 -o ${INPUT_VCF%.*.*}.comp.vcf.gz $INPUT_VCF \
    && rm ${INPUT_VCF}*
