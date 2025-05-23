#!/bin/bash
# Set up to run samples or workfile
# First argument is either a directory named with the FULLSMID that holds FQs, or
# a workfile of directories with the last directory being named with the FULLSMID
if [[ -f $1 ]]; then INDIRS=($(cat $1)); else INDIRS=($1); fi
export THREADS=8
# Second argument is the directory where all samples will end up in
RUN_OUTDIR=$2
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
for INDIR in ${INDIRS[@]}; do

##0. Set up environment for processing
FULLSMID=${INDIR##*/}
ENV_FILE=${INDIR}/${FULLSMID}.env
echo -e "ENV_FILE=${ENV_FILE}" > $ENV_FILE
echo -e "FULLSMID=${FULLSMID}" >> $ENV_FILE
echo -e "INDIR=${INDIR}" >> $ENV_FILE
[ ! -d ${RUN_OUTDIR} ] && mkdir ${RUN_OUTDIR}
echo -e "OUTDIR=${RUN_OUTDIR}/${FULLSMID}" >> $ENV_FILE
[ ! -d ${RUN_OUTDIR}/${FULLSMID} ] && mkdir ${RUN_OUTDIR}/${FULLSMID}
echo -e "LOG_FILE=${RUN_OUTDIR}/${FULLSMID}/${FULLSMID}.log" >> $ENV_FILE
echo -e "RUN_TYPE=exome" >> $ENV_FILE
echo -e "BAM=${FULLSMID}.aln.srt.mrk.bam" >> $ENV_FILE
echo -e "CRAM=${FULLSMID}.aln.srt.mrk.cram" >> $ENV_FILE
echo -e "GVCF=${FULLSMID}.snp.indel.g.vcf.gz" >> $ENV_FILE
echo -e "STATS_FILE=${RUN_OUTDIR}/${FULLSMID}/${FULLSMID}.stats.csv" >> $ENV_FILE
echo -e "TMP_DIR=${HOME}/WEStmp" >> $ENV_FILE
[ ! -d ${HOME}/WEStmp ] && mkdir ${HOME}/WEStmp
cat ${SCRIPT_DIR}/../baseEnvs/references_hydra.env >> $ENV_FILE

for VAR in $(cat $ENV_FILE); do export $VAR; done

cp $ENV_FILE $OUTDIR
echo -n "" > $LOG_FILE

## 1. Align and Sort
for FQ in $(find $INDIR -name "*1.f*q.gz"); do
SM=$(echo $FULLSMID | cut -d^ -f1)
BARCODE=$(echo $FULLSMID | cut -d^ -f2)
PROJECT=$(echo $FULLSMID | cut -d^ -f3)
HEAD=$(zcat ${FQ} | head -n1)
HEAD_CHECK="${HEAD//[^:]}"
if [[ ${#HEAD_CHECK} -eq 4 ]]; then
FLOWCELL=$(zcat ${FQ} | head -n1 | cut -d: -f1 | cut -d '@' -f2)
LANE=$(zcat ${FQ} | head -n1 | cut -d: -f2)
FLOWLANE="${FLOWCELL}.${LANE}"
elif [[ ${#HEAD_CHECK} -gt 4 ]]; then
FLOWCELL=$(zcat ${FQ} | head -n1 | cut -d: -f3)
LANE=$(zcat ${FQ} | head -n1 | cut -d: -f4)
FLOWLANE="${FLOWCELL}.${LANE}"
else
FLOWLANE=$(echo ${FQ##*/} | rev | cut -d_ -f2- | rev)
fi

echo "@RG\tID:${FLOWLANE}\tPL:illumina\tPU:${FLOWLANE}:${BARCODE}\tLB:${BARCODE}\tSM:${SM}\tDS:${FULLSMID}" > ${OUTDIR}/${FULLSMID}.${FLOWLANE}.rgfile

RG="${OUTDIR}/${FULLSMID}.${FLOWLANE}.rgfile"
bwa-mem2 mem -M -t $THREADS -K 10000000 \
  -R $(head -n1 $RG) \
  ${REF_FASTA} \
  ${FQ} \
  ${FQ/1.f/2.f} \
  | gatk \
  --java-options "-Xmx70g -XX:ParallelGCThreads=1" \
  SortSam \
  -I /dev/stdin \
  -O ${OUTDIR}/${FQ##*/}.bam \
  -R ${REF_FASTA} \
  -SO coordinate \
  --MAX_RECORDS_IN_RAM 1000000 \
  --CREATE_INDEX true \
  2>> $LOG_FILE

# 1.2 Extract exome by intersecting the aligned bam
bedtools intersect -u -a ${OUTDIR}/${FQ##*/}.bam -b $REF_PADBED > ${OUTDIR}/${FQ##*/}.isec.bam \
&& rm ${OUTDIR}/${FQ##*/}.ba* \
2>> $LOG_FILE

done

# 1.3 Mark Duplicates
MD_INPUTS=()
for BM in $(find $OUTDIR -name "*.isec.bam"); do
MD_INPUTS+="-I ${BM} "
done
gatk \
  --java-options "-Xmx40g -XX:ParallelGCThreads=1" \
  MarkDuplicates \
    ${MD_INPUTS[@]}\
    -O ${OUTDIR}/${BAM} \
    -M ${OUTDIR}/${FULLSMID}.dup.metrics.txt \
    -R ${REF_FASTA} \
    --TMP_DIR ${TMP_DIR} \
    --CREATE_INDEX true \
    2>> $LOG_FILE
#rm -R ${INDIR}
rm $OUTDIR/*isec.ba*

## 2. BQSR - Recalibrate Bases
# 2.1 Generate Recal Table
${GATK} \
  --java-options "-Xmx40g -XX:ParallelGCThreads=1" \
  BaseRecalibratorSpark \
    -I ${OUTDIR}/${BAM} \
    -R ${REF_FASTA} \
    -L ${REF_PADBED} \
    --known-sites ${REF_MILLS_GOLD} \
    --known-sites ${REF_DBSNP} \
    --known-sites ${REF_ONEKGP1} \
    -O ${OUTDIR}/${FULLSMID}.recal.txt \
    -- \
    --spark-master local[$THREADS] \
    2>> $LOG_FILE

# 2.2 Apply Recal Table
gatk \
  --java-options "-Xmx40g -XX:ParallelGCThreads=1" \
  ApplyBQSR \
    -I ${OUTDIR}/${BAM} \
    -bqsr ${OUTDIR}/${FULLSMID}.recal.txt \
    -R ${REF_FASTA} \
    -L ${REF_PADBED} \
    -O ${OUTDIR}/${FULLSMID}.recal.bam \
    2>> $LOG_FILE

samtools index ${OUTDIR}/${FULLSMID}.recal.bam

## 3. Call Variants
gatk \
  --java-options "-Xmx40g -XX:ParallelGCThreads=1" \
  HaplotypeCaller \
    -I ${OUTDIR}/${FULLSMID}.recal.bam \
    -R ${REF_FASTA} \
    --dbsnp ${REF_DBSNP} \
    -L ${REF_PADBED} \
    -ERC GVCF \
    -O ${OUTDIR}/${GVCF} \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    2>> $LOG_FILE

rm ${OUTDIR}/${FULLSMID}.recal.ba*

## 5. QC
#5.1 Coverage
#5.1a Raw WES Coverage
gatk \
  --java-options "-Xmx20g -XX:ParallelGCThreads=1" \
  DepthOfCoverage \
    -I ${OUTDIR}/${BAM} \
    -O ${OUTDIR}/${BAM}.rawwgsmetrics.txt \
    -L ${REF_PADBED%.bed}.interval_list \
    -R ${REF_FASTA} \
    --summary-coverage-threshold 10 \
    --summary-coverage-threshold 15 \
    --summary-coverage-threshold 20 \
    --summary-coverage-threshold 30 \
    --summary-coverage-threshold 40 \
    --summary-coverage-threshold 50 \
    --summary-coverage-threshold 60 \
    --summary-coverage-threshold 70 \
    --summary-coverage-threshold 80 \
    --summary-coverage-threshold 90 \
    --summary-coverage-threshold 100 \
    --omit-depth-output-at-each-base true \
    --omit-interval-statistics true \
    --omit-locus-table true

#5.2 FREEMIX
VerifyBamID2 \
  --BamFile ${OUTDIR}/${BAM} \
  --SVDPrefix ${REF_VBID} \
  --Reference ${REF_FASTA} \
  --NumThread ${THREADS} \
  --Output ${OUTDIR}/${BAM}.vbid2 \
  --max-depth 1000

#5.3 Variant Calling Metrics
gatk \
  --java-options "-Xmx30g -XX:ParallelGCThreads=1" \
  CollectVariantCallingMetrics \
    -I ${OUTDIR}/${GVCF} \
    -O ${OUTDIR}/${GVCF##*/}.vcfmetrics \
    -TI ${REF_PADBED%.bed}.interval_list \
    -R ${REF_FASTA} \
    --DBSNP ${REF_DBSNP} \
    --THREAD_COUNT 6 \
    --GVCF_INPUT true \
    --TMP_DIR ${TMP_DIR}

#CREATE CRAM
samtools view -C -T $REF_FASTA --output-fmt-option archive --output-fmt-option embed_ref -o $OUTDIR/$CRAM $OUTDIR/$BAM \
&& samtools index $OUTDIR/$CRAM \
&& rm $OUTDIR/$BAM

#5.4 SNPeff Annotations
# Creates Fields file local to gVCF
FIELDS_FILE=${OUTDIR}/${GVCF}.snpeff-5.1-FIELDS.txt

## 1st annotate with SNPEFF
bcftools view -R ${KEYGENES} ${OUTDIR}/${GVCF} \
  | java -Xmx20g -jar ${SNPEFF} ann -interval ${KEYGENES} -noStats -v GRCh38.105 \
  | /scripts/snpEff_v5.1/scripts/vcfEffOnePerLine.pl \
  | java -Xmx20g -jar ${SNPSIFT} extractFields -e "."  \
  - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" \
  "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" \
  > ${FIELDS_FILE}

  [ -z $STATS_FILE ] && STATS_FILE=${OUTDIR}/${FULLSMID}.stats.csv

  ## RG comparison
  IN_RGS=$(ls ${OUTDIR} | grep -c rgfile)
  CRAM_RGS=$(samtools view -H ${OUTDIR}/${FULLSMID}*.cram | grep -wc '^@RG')

  ## Get Coverage Stats
  RAW_COV=$(head -n8 ${OUTDIR}/${FULLSMID}*.cram.rawwgsmetrics.txt | tail -n1 | cut -f2)
  WGS_COV=$(head -n8 ${OUTDIR}/${FULLSMID}*.cram.wgsmetrics.txt | tail -n1 | cut -f2)
  PADEX_COV=$(head -n8 ${OUTDIR}/${FULLSMID}*.cram.wgsmetrics_paddedexome.txt | tail -n1 | cut -f2)

  ## Get FREEMIX
  FREEMIX=$(tail -n1 ${OUTDIR}/${FULLSMID}*.cram.vbid2.selfSM | cut -f7)

  ## Get TITV Ratio
  TOTAL_SNPS=$(head -n8 ${OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f1)
  TOTAL_INDELS=$(head -n8 ${OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f8)
  PCT_DBSNP=$(head -n8 ${OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f5)
  NOVEL_TITV=$(head -n8 ${OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f7)
  DBSNP_TITV=$(head -n8 ${OUTDIR}/${FULLSMID}*.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics | tail -n1 | cut -f6)
  TOTAL_TITV=$(bc -l <<<"($PCT_DBSNP*$DBSNP_TITV)+((1-$PCT_DBSNP)*$NOVEL_TITV)")

  ## Get Annotation stats
  FIELDS_FILE="${OUTDIR}/${FULLSMID}*.vcf.gz.snpeff-5.1-FIELDS.txt"
  GENES=(APP PSEN1 PSEN2 GRN MAPT MAPT MAPT TREM2)
  CAN_TS=(ENST00000346798.8 ENST00000324501 ENST00000366783.8	 ENST00000053867.8 ENST00000262410.10 ENST00000629368.2 ENST00000618029.3 ENST00000373113.8)
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

done
