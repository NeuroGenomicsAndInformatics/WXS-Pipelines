#!/bin/bash
# This script adds standard annotations to a vcf file.
INPUT_DIR=$1
INTERVAL=$2
[[ -z $INTERVAL ]] && INTERVAL=${LSB_JOBINDEX}
INPUT_VCF=$(find ${INPUT_DIR} -name "*.${INTERVAL}.*.splitMA.geno.sites.vcf.gz")

INT_BEDS=($(ls ${INPUT_DIR}/intlists | grep scattered.bed))
INT_BED="${INPUT_DIR}/intlists/${INT_BEDS[${INTERVAL}-1]}"
CHR=$(echo ${INPUT_DIR##*/} | cut -d_ -f1)

SNPEFF="/ref/snpEff_v5.1/snpEff.jar"
SNPSIFT="/ref/snpEff_v5.1/SnpSift.jar"

OUTPUT_VCF="${INPUT_VCF%.*.*}.ann.vcf.gz"

java -Xmx20g -jar ${SNPEFF} \
	ann \
	-canon \
	-interval ${INT_BED} \
	-noStats \
	-v GRCh38.105 \
	$INPUT_VCF \
	| bgzip -c > ${INPUT_VCF%.*.*}.snpeff.vcf.gz

java -Xmx20g -jar ${SNPSIFT} annotate \
	-id \
	-info dbSNPBuildID \
	/ref/ADSPR5/dbSNP_b156_GRCh38.gz \
	-v ${INPUT_VCF%.*.*}.snpeff.vcf.gz \
	| bgzip -c > ${INPUT_VCF%.*.*}.snpeff.dbSNP.vcf.gz

java -Xmx20g -jar ${SNPSIFT} annotate \
	-id \
	-info AF,AC,AF_popmax,AC_popmax,AF_nfe,AC_nfe,AF_afr,AC_afr,AF_amr,AC_amr,AF_sas,AC_sas \
	/ref/gnomad/gnomad.genomes.v4.1.sites.${CHR}.vcf.bgz \
	-v ${INPUT_VCF%.*.*}.snpeff.dbSNP.vcf.gz \
	| bgzip -c > ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.vcf.gz
	
bcftools annotate \
	--rename-chrs /ref/ADSPR5/chr_change.txt \
	-Oz \
	-o ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.vcf.gz \
	${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.vcf.gz
	
bcftools annotate \
	-a /ref/ADSPR5/whole_genome_SNVs.tsv.gz \
	-c CHROM,POS,REF,ALT,CADD_RAW,CADD_PHRED \
	-h /ref/ADSPR5/CADD_hdr.txt \
	-Oz \
	-o ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.vcf.gz \
	${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.vcf.gz	

bcftools annotate \
	-a /ref/ADSPR5/revel_with_transcript_ids_fixed.tsv.gz \
	-c CHROM,POS,REF,ALT,REVEL_SCORE \
	-h /ref/ADSPR5/revel_hdr.txt \
	-Oz \
	-o ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.vcf.gz \
	${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.vcf.gz	
tabix -f ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.vcf.gz

bcftools annotate \
	-a /ref/ADSPR5/clinvar_20240708.vcf.gz \
	-c CLNSIG \
	${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.vcf.gz \
	| bgzip -c > ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.clinvar.vcf.gz
	
bcftools annotate \
	--rename-chrs /ref/ADSPR5/chr_changeback.txt \
	-Oz \
	-o ${OUTPUT_VCF} \
	${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.clinvar.vcf.gz
	
SUCCESS=$?
[[ $SUCCESS -eq 0 ]] && [[ -s ${OUTPUT_VCF} ]] \
	&& rm ${INPUT_VCF%.*.*}.snpeff.vcf.gz \
	&& rm ${INPUT_VCF%.*.*}.snpeff.dbSNP.vcf.gz \
	&& rm ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.vcf.gz \
	&& rm ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.vcf.gz \
	&& rm ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.vcf.gz \
	&& rm ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.vcf.gz \
	&& rm ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.vcf.gz.tbi \
	&& rm ${INPUT_VCF%.*.*}.snpeff.dbSNP.gnomad.chrs.CADD.REVEL.clinvar.vcf.gz \
	&& rm ${INPUT_VCF} \
	&& rm ${INPUT_VCF}.tbi \

exit $SUCCESS