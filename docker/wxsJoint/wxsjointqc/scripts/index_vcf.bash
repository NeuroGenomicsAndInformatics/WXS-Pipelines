#!/bin/bash
# This script adds a missing index to a vcf file.
INPUT_VCF=$1

[[ -f ${INPUT_VCF}.tbi ]] || tabix ${INPUT_VCF}
