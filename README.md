This WXS-Pipelines repository is designed for processing whole genome
and exome sequencing data on WUSTL's compute1 HPC cluster.

Inputs:
1. A text file with sample names (FULLSMID) separated by whitespace. Spaces, tabs, or newlines are acceptable.
  ex. MAP_10040.GRAYS^2102557103^20202110_MGI_WGS_GRAYS
2. Data staged in storage1 with the following structure:
  /storage1/fs1/cruchagac/Active/$USER/c1in/${FULLSMID}/${FLOWCELL}.${LANE}_${READ}.fastq.gz
  ex. /storage1/fs1/cruchagac/Active/matthewj/c1in/MAP_10040^2102557103^20202110_MGI_WGS_GRAYS/HLFK3DSXY.1_1.fastq.gz
  /storage1/fs1/cruchagac/Active/matthewj/c1in/MAP_10040^2102557103^20202110_MGI_WGS_GRAYS/HLFK3DSXY.1_2.fastq.gz
  /storage1/fs1/cruchagac/Active/matthewj/c1in/MAP_10040^2102557103^20202110_MGI_WGS_GRAYS/HLNHFDSXY.3_1.fastq.gz
  /storage1/fs1/cruchagac/Active/matthewj/c1in/MAP_10040^2102557103^20202110_MGI_WGS_GRAYS/HLNHFDSXY.3_2.fastq.gz
  ...
3. Fastqs need specific filenames with the flowcell and lane for RGfile generation, but bams and crams can be named however as long as they're in the right directory.

Outputs:
1. Aligned, sorted, duplicate-marked crams for each sample
  ex. MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.aln.srt.mrk.cram
      MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.aln.srt.mrk.cram.crai
2. A gVCF for the sample
  ex. MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.snp.indel.g.vcf.gz
      MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.snp.indel.g.vcf.gz.tbi
3. Coverage reports
  ex. MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.aln.srt.mrk.cram.rawwgsmetrics.txt
      MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.aln.srt.mrk.cram.wgsmetrics.txt
      MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.aln.srt.mrk.cram.wgsmetrics_paddedexome.txt
4. VerifyBamID report
  ex. MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.aln.srt.mrk.cram.vbid2.Ancestry
      MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.aln.srt.mrk.cram.vbid2.selfSM
5. Variant Calling Metrics
  ex. MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.snp.indel.g.vcf.gz.vcfmetrics.variant_calling_detail_metrics
      MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.snp.indel.g.vcf.gz.vcfmetrics.variant_calling_summary_metrics
6. Various outputs, logs, and reports
   ex. MAP_10040^2102557103^20202110_MGI_WGS_GRAYS.env - The environment variables used to process the sample


STEPS TO RUN PIPELINE
1. Stage data to storage1
  a. Use the structure defined above in "Inputs"
  b. There is also a script in the repository called prepareDataForStaging.bash that demonstrates how one might stage with a script and can be modified.
2. Log into compute1 (if you aren't already)
    ssh USER@compute1-client-1.ris.wustl.edu
    NOTE: The above uses client 1, but there are 4 client nodes numbered 1-4
3. Create a workfile
  a. Any text file should work. It should contain a list of FULLSMIDs that correspond to the folder names of the staged data.
  b. Spaces, tabs, or new lines for each file should work and have all been tested.
  c. Running a list command to get the folders is a quick way of doing this.
    ls /storage1/fs1/cruchagac/matthewj/c1in/*GRAYS > workfile.txt
4. Get the pipeline code
  a. The pipelines are publicly available via github. You can get them by navigating to where you want the folder and entering:
    git clone https://github.com/NeuroGenomicsAndInformatics/WXS-Pipelines.git
5. Run the pipeline
  a. Navigate to the folder that was just added.
    cd /path/to/WXS-Pipelines
  b. Call the pipeline script
    ex. bash wgsfqspipeline.bash /path/to/workfile.txt
6. Collect outputs when done
  a. The output files will show up in /storage1/fs1/cruchagac/$USER/c1out/$FULLSMID
  b. The pipeline repository also includes a transferOutputFiles.bash script to copy files back to storage1. The input could also be the same workfile.

The steps in the processing pipeline are below:
1. Prepare fastq files (if necessary)
  a. Crams are first reverted with RevertSam to remove all prior tags to ensure clean processing from scratch. The reverted data is piped to SamToFastq.
  b. SamToFastq separates the data into paired-end fastq files by read group. This step ensures that all samples are processed identically by creating a common state for the raw data. This common state is read group-separated fastq files.
2. Align and sort fastq reads
  a. Reference files for the alignments, and the entire pipeline, were taken from GATK’s resource bundle and were downloaded from the Google bucket on August 12 of 2019 and are GRCh38. 
  b. The reads were aligned with either bwa-mem2 on CPUs or NVIDIA’s Clara Parabricks version 3.6.1 on GPUs. Both of these perform the same function and produce the same outputs as bwa-mem version 0.7.15 and 0.7.17. The alignments were done in alternate-aware mode with alternate contigs in the GATK resource bundle. Bwa-mem2 also requires index files that were locally generated from the resource bundle reference fasta.
  c. The reads are coordinate-sorted with GATK’s SortSam tool. This tool is internally a Picard tool.
3. Mark duplicates
  a. The marking of duplicate reads was done with either GATK 4.2.6.1 or Parabricks 3.6.1. These tools perform the same functions and lead to the same outcomes, but the Parabricks version is able to take advantage of GPUs to massively parallelize the job. In the case that a GPU wasn’t available, GATK performed the work on CPUs.
  b. We are able to utilize the read data in cram format beyond this point. After marking duplicates, the read data is stored in cram format and those are used to process the data further.
4. Recalibrate bases
  a. GATK BaseRecalibratorSpark was utilized to perform base recalibration in parallel. While providing identical results, this tool was able to recalibrate base and mapping qualities an order of magnitude faster than its serial counterpart.
  b. Mills-Gold
  c. dbSNP 138
  d. 1K genome Phase 1
5. Call variants with HaplotypeCaller
  a. Parabricks 3.6.1 was utilized for variant calling on the genomes. Parabricks 3.6.1 performs the HaplotypeCaller algorithm from GATK 4.1.0 with near identical results. It also saves the step of applying the BQSR values to a bam. Because of the scale of data being processed, Parabricks’ implementation was the preferred option because others were prohibitively slow.
QC components:
1. Mean coverage
  a. GATK 4.2.6.1’s CollectRawWgsMetrics and CollectWgsMetrics provided coverage data for the samples. Summary data was generated for coverage at base and mapping quality of 3 and 0, respectively, by CollectRawWgsMetrics. CollectWgsMetrics represents a higher quality threshold for the data with base and mapping quality minimums of 20, corresponding to 99% certainty.
  b. Coverage data was also collected with CollectWgsMetrics over only the padded exome.
2. Freemix
  a. Freemix was calculated with Griffan Lab’s VerifyBamID, an updated version of the original. These are used to measure the probability of contamination.
  b. The data used for this analysis come from the 1000 Genomes – Phase 3 data release. The files used can be found in the repository for the tool.
3. TiTv Ratio
  a. GATK 4.2.6.1’s CollectVariantCallingMetrics was used to capture a variety of data from the gVCFs. Of note are the Percent SNPS in dbSNP, dbSNP TiTv Ratio, and Novel TiTv Ratio that we use to measure the total TiTv ratio for the data.
4. Annotations of key genes
  a. SnpEff/SnpSift 5.1 annotate APP, PSEN1, PSEN2, GRN, MAPT, TREM2
5. Stats File
  a. A csv with a header line and a data line.
  b. Key QC stats from above
  c. Counts of annotations on MANE transcripts from the key genes above
  d. Counts of HIGH and MODERATE effect annotations from the key genes above
  e. The HIGH and MODERATE effect variants from the key genes above  