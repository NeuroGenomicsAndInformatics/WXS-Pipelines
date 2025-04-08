#!/usr/bin/env -S Rscript --vanilla
# Script to collect counts files from intervals

library(data.table)

args <- commandArgs(TRUE)
input_dir <- args[1]
namebase = basename(input_dir)

count_files <- list.files(path = input_dir, pattern = "*.counts.csv", full.names = TRUE)
combined_counts <- ""
for (count_file in count_files) {
  if (count_file == count_files[1]) {
    combined_counts <- read.csv(count_file)
  } else {
    count_csv <- read.csv(count_file)
    combined_counts <- rbind(combined_counts,count_csv)
  }	
}

combined_counts[nrow(combined_counts)+1,] <- c("Totals",sum(combined_counts$Original),sum(combined_counts$QUAL),sum(combined_counts$LCR),sum(combined_counts$no_mono),sum(combined_counts$maxDP),sum(combined_counts$ABHet),sum(combined_counts$splitMA),sum(combined_counts$geno))

write.csv(combined_counts, file = paste(input_dir,"/",namebase,"_combined_counts.csv", sep=""), row.names = FALSE)