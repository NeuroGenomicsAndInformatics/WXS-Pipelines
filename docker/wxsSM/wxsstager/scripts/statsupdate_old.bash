#!/bin/bash
[ -z $FULLSMID ] && FULLSMID=${1##*/}
[ -z $FINAL_OUTDIR ] && FINAL_OUTDIR=$1
[ -z $STATS_FILE ] && STATS_FILE=${FINAL_OUTDIR}/${FULLSMID}.stats.csv

## RG comparison
IN_RGS=$(ls ${FINAL_OUTDIR} | grep -c rgfile)
CRAM_RGS=$(samtools view -H ${FINAL_OUTDIR}/${FULLSMID}*.cram | grep -wc '^@RG')

## Get Coverage Stats
RAW_COV=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.cram.rawwgsmetrics.txt | tail -n1 | cut -f2)
WGS_COV=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.cram.wgsmetrics.txt | tail -n1 | cut -f2)
PADEX_COV=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.cram.wgsmetrics_paddedexome.txt | tail -n1 | cut -f2)

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
GENES=(APP PSEN1 PSEN2 GRN MAPT MAPT MAPT TREM2)
CAN_TS=(ENST00000346798.8 ENST00000324501 ENST00000366783.8 ENST00000053867.8 ENST00000262410.10 ENST00000629368.2 ENST00000618029.3 ENST00000373113.8)
NUM_ANNS=()
NUM_ANNS_HM=()
for TS in ${CAN_TS[@]}; do
  NUM_ANNS+=($(grep -c ${TS} ${FIELDS_FILE}))
  NUM_ANNS_HM+=($(grep ${TS} ${FIELDS_FILE} | grep -c -E 'HIGH|MODERATE'))
done

## Set Header Line
echo -n "FULLSMID,RGs in,Cram RGs,Raw Mean Coverage,Mean Coverage,Padded Exome Mean Coverage,FREEMIX,Total TITV,PCT dbSNP,dbSNP TITV,Novel TITV,Num SNPs,Num Indels," > ${STATS_FILE}
for GENE in ${GENES[@]}; do echo -n "$GENE Mane Annotations,$GENE HM," >> ${STATS_FILE}; done
echo "" >> ${STATS_FILE}
## Add stats
echo -n "${FULLSMID},${IN_RGS},${CRAM_RGS},${RAW_COV},${WGS_COV},${PADEX_COV},${FREEMIX},${TOTAL_TITV},${PCT_DBSNP},${DBSNP_TITV},${NOVEL_TITV},${TOTAL_SNPS},${TOTAL_INDELS}," >> ${STATS_FILE}
for ((i = 0 ; i < ${#CAN_TS[@]} ; i++)); do
echo -n "${NUM_ANNS[$i]},${NUM_ANNS_HM[$i]}," >> ${STATS_FILE}; done
echo "" >> ${STATS_FILE}
