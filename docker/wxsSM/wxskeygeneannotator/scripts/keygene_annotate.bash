#!/bin/bash
#To run below "pipeline" in roughly 5 samples per sequencing project to have an idea of whether all projects used the wrong annotation files or only some
# First argument is a text file with absolute or relative paths to gVCFs to check or a single gVCF
GVCF=$1

# Second optional argument is for an output CSV to store stats. Default is output.csv in current directory
OUTDIR=${GVCF%/*}
FULLSMID=$(echo "${GVCF##*/}" | cut -d '.' -f1)
OUTFILE=${OUTDIR}/${FULLSMID}.keygenes.csv

# Set locations for SNPEFF and SNPSIFT
SNPEFF="/scripts/snpEff_v5.1/snpEff.jar"
SNPSIFT="/scripts/snpEff_v5.1/SnpSift.jar"
KEYGENES="/scripts/keygenes.bed"

# Create output file with header
echo "FULLSMID,MAPT Annotations,MAPT Annotations in MANE,MAPT MANE HIGH OR MODERATE\
,APP Annotations,APP Annotations in MANE,APP MANE HIGH OR MODERATE\
,PSEN1 Annotations,PSEN1 Annotations in MANE,PSEN1 MANE HIGH OR MODERATE\
,PSEN2 Annotations,PSEN2 Annotations in MANE,PSEN2 MANE HIGH OR MODERATE" > ${OUTFILE}


# Run MAPT check
# Supports single file or workfile
# Creates Fields file local to gVCF
FIELDS_FILE=${GVCF%.*}-snpeff-5.1-FIELDS.txt

## 1st annotate with SNPEFF
bcftools view -R ${KEYGENES} ${GVCF} \
  | java -Xmx20g -jar ${SNPEFF} ann -interval ${KEYGENES} -noStats -v GRCh38.105 \
  | /scripts/snpEff_v5.1/scripts/vcfEffOnePerLine.pl \
  | java -Xmx20g -jar ${SNPSIFT} extractFields -e "."  \
  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \
  "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" \
  > ${FIELDS_FILE}
  
# Add to output file
echo "$FULLSMID\
,$(grep -wc MAPT ${FIELDS_FILE}),$(grep -wc ENST00000262410.10 ${FIELDS_FILE}),$(grep ENST00000262410.10 ${FIELDS_FILE} | grep -c -E "HIGH|MODERATE") \
,$(grep -wc APP ${FIELDS_FILE}),$(grep -wc ENST00000346798.8 ${FIELDS_FILE}),$(grep ENST00000346798.8 ${FIELDS_FILE} | grep -c -E "HIGH|MODERATE") \
,$(grep -wc PSEN1 ${FIELDS_FILE}),$(grep -wc ENST00000324501.10 ${FIELDS_FILE}),$(grep ENST00000324501.10 ${FIELDS_FILE} | grep -c -E "HIGH|MODERATE") \
,$(grep -wc PSEN2 ${FIELDS_FILE}),$(grep -wc ENST00000366783.8 ${FIELDS_FILE}),$(grep ENST00000366783.8 ${FIELDS_FILE} | grep -c -E "HIGH|MODERATE")" \
>> $OUTFILE
