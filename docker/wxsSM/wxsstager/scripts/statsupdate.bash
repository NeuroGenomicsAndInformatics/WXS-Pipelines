#!/bin/bash
[ -z $FULLSMID ] && FULLSMID=$1
[ -z $STAGE_INDIR ] && STAGE_INDIR=/storage1/fs1/cruchagac/Active/$USER/c1in/$FULLSMID
[ -z $OUTDIR ] && OUTDIR=/scratch1/fs1/cruchagac/$USER/c1out/$FULLSMID
[ -z $FINAL_OUTDIR ] && FINAL_OUTDIR=/storage1/fs1/cruchagac/Active/$USER/c1out/$FULLSMID
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
NOVEL_TITV=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f6)
DBSNP_TITV=$(head -n8 ${FINAL_OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f7)
TOTAL_TITV=$(( (PCT_DBSNP * DBSNP_TITV) + ((1-PCT_DBSNP) * NOVEL_TITV) ))

## Get Annotation stats
FIELDS_FILE="${FINAL_OUTDIR}/${FULLSMID}*.vcf-snpeff-5.1-FIELDS.txt"
GENES=(APOE APP ATP13A2 CHMP2B DNAJC16 EIF4G1 FBXO7 FUS GIGYF2 GRN HTRA2 LRRK2 MAPT MAPT MAPT OPTN PARK7 PFN1 PINK1 PLA2G6 PSEN1 PSEN2 SNCA SOD1 SQSTM1 SQSTM1 TARDBP TBK1 UBQLN2 UCHL1 VCP VPS35)
CAN_TS=(ENST00000252486.9 ENST00000346798.8 ENST00000452699.5 ENST00000263780.9 ENST00000375847.8 ENST00000342981.8 ENST00000266087.12 ENST00000254108.12 ENST00000373563.9 ENST00000053867.8 ENST00000258080.8 ENST00000298910.12 ENST00000262410.10 ENST00000626571.2 ENST00000620818.4 ENST00000378747.8 ENST00000493678.5 ENST00000225655.6 ENST00000321556.5 ENST00000660610.1 ENST00000700267.1 ENST00000422240.6 ENST00000336904.7 ENST00000270142.11 ENST00000389805.9 ENST00000643389.2 ENST00000240185.8 ENST00000331710.10 ENST00000338222.7 ENST00000284440.9 ENST00000358901.11 ENST00000299138.12)
NUM_ANNS=()
NUM_ANNS_HM=()
for TS in ${CAN_TS[@]}; do
  NUM_ANNS+=($(grep -c ${TS} ${FIELDS_FILE}))
  NUM_ANNS_HM+=($(grep ${TS} ${FIELDS_FILE} | grep -c -E 'HIGH|MODERATE'))
done

## Set Header Line
echo -n "FULLSMID,RGs in,Cram RGs,Raw Mean Coverage,Mean Coverage,Padded Exome Mean Coverage,FREEMIX,Total TITV,dbSNP TITV,Novel TITV,Num SNPs,Num Indels," > ${STATS_FILE}
for GENE in ${GENES[@]}; do echo -n "$GENE Mane Annotations,$GENE HM," >> ${STATS_FILE}; done
echo "" >> ${STATS_FILE}
## Add stats
echo -n "${FULLSMID},${IN_RGS},${CRAM_RGS},${RAW_COV},${WGS_COV},${PADEX_COV},${FREEMIX},${TOTAL_TITV},${DBSNP_TITV},${NOVEL_TITV},${TOTAL_SNPS},${TOTAL_INDELS}," >> ${STATS_FILE}
for ((i = 0 ; i < ${#CAN_TS[@]} ; i++)); do
echo -n "${NUM_ANNS[$i]},${NUM_ANNS_HM[$i]}," >> ${STATS_FILE}; done
echo "" >> ${STATS_FILE}
