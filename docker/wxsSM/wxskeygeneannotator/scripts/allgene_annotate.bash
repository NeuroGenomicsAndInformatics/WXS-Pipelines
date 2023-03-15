#!/bin/bash
#To run below "pipeline" in roughly 5 samples per sequencing project to have an idea of whether all projects used the wrong annotation files or only some
# First argument is a text file with absolute or relative paths to gVCFs to check or a single gVCF
if [[ -z $FULLSMID ]]; then FULL_GVCF=$1; FINAL_OUTDIR=${FULL_GVCF%/*}; GVCF=${FULL_GVCF##*/}; FULLSMID=$(echo "${GVCF}" | cut -d '.' -f1); fi

# Set locations for SNPEFF and SNPSIFT
SNPEFF="/scripts/snpEff_v5.1/snpEff.jar"
SNPSIFT="/scripts/snpEff_v5.1/SnpSift.jar"

# Creates Fields file local to gVCF
FIELDS_FILE=${FINAL_OUTDIR}/${GVCF}.snpeff-5.1-FIELDS.txt

## 1st annotate with SNPEFF
java -Xmx20g -jar ${SNPEFF} ann -noStats -v GRCh38.105 ${FINAL_OUTDIR}/${GVCF} \
  | /scripts/snpEff_v5.1/scripts/vcfEffOnePerLine.pl \
  | java -Xmx20g -jar ${SNPSIFT} extractFields -e "."  \
  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \
  "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" \
  > ${FIELDS_FILE}
