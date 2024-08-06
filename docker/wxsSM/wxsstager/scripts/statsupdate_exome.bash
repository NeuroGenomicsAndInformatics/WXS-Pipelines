#!/bin/bash
[ -z $FULLSMID ] && FULLSMID=${1##*/}
[ -z $FINAL_OUTDIR ] && FINAL_OUTDIR=$1
[ -z $STATS_FILE ] && STATS_FILE=${FINAL_OUTDIR}/${FULLSMID}.stats.csv

## RG comparison
IN_RGS=$(ls ${FINAL_OUTDIR} | grep -c rgfile)
CRAM_RGS=$(samtools view -H ${FINAL_OUTDIR}/${FULLSMID}*.cram | grep -wc '^@RG')

## Get Coverage Stats
PADEX_COV=$(head -n2 ${FINAL_OUTDIR}/${FULLSMID}*.cram.docmetrics_paddedexome.sample_summary | tail -n1 | cut -d, -f3)

## Get FREEMIX
FREEMIX=$(tail -n1 ${FINAL_OUTDIR}/${FULLSMID}*.cram.vbid2.selfSM | cut -f7)

## Get TITV Ratio
TOTAL_SNPS=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f1)
TOTAL_INDELS=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f8)
PCT_DBSNP=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f5)
NOVEL_TITV=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f7)
DBSNP_TITV=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f6)
TOTAL_TITV=$(bc -l <<<"($PCT_DBSNP*$DBSNP_TITV)+((1-$PCT_DBSNP)*$NOVEL_TITV)")

## Get Annotation stats
FIELDS_FILE="${FINAL_OUTDIR}/${FULLSMID}*.vcf.gz.snpeff-5.1-FIELDS.txt"
GENES=(APP PSEN1 PSEN2 GRN MAPT TREM2 SNCA LRRK2 PARK2 PINK1)
CAN_TS=(ENST00000346798.8 ENST00000324501 ENST00000366783.8 ENST00000053867.8 ENST00000262410.10 ENST00000373113.8 ENST00000394991 ENST00000298910 ENST00000366898 ENST00000321556)

## Set Header Line
echo -n "FULLSMID,RGs in,Cram RGs,Padded Exome Mean Coverage,FREEMIX,Total TITV,PCT dbSNP,dbSNP TITV,Novel TITV,Num SNPs,Num Indels," > ${STATS_FILE}
for GENE in ${GENES[@]}; do echo -n "$GENE Mane Annotations count,$GENE HM count,$GENE Annotations," >> ${STATS_FILE}; done
echo "" >> ${STATS_FILE}
## Add stats
echo -n "${FULLSMID},${IN_RGS},${CRAM_RGS},${PADEX_COV},${FREEMIX},${TOTAL_TITV},${PCT_DBSNP},${DBSNP_TITV},${NOVEL_TITV},${TOTAL_SNPS},${TOTAL_INDELS}," >> ${STATS_FILE}
## Add Annotations
for TS in ${CAN_TS[@]}; do
echo -n "$(grep -c ${TS} ${FIELDS_FILE}),$(grep ${TS} ${FIELDS_FILE} | grep -c -E 'HIGH|MODERATE')," >> ${STATS_FILE}
for VAR in $(grep ${TS} ${FIELDS_FILE} | grep -E 'HIGH|MODERATE' | cut -f14); do echo -n "$VAR " >> ${STATS_FILE}; done
echo -n "," >> ${STATS_FILE}
done
echo "" >> ${STATS_FILE}

## Get sex information
bash /scripts/sex_check_coverage.bash ${FINAL_OUTDIR}/${FULLSMID}*.cram