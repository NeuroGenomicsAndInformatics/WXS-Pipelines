#!/usr/bin/env -S Rscript --vanilla --slave

### 2017-03-01 - vifehe 
### Small script to calcualte average DP from  GATK extracted table
### UAGE: ./CALCULATE_DP.R ${INFILE}


library(tools)
require(data.table)
library(data.table)

cat("\nreading passed by arguments\n")
args <- commandArgs(TRUE)
INFILE<-args[1]


#Open R
library(tools)
tables<-Sys.glob(paste0(INFILE))
for (TABLE in 1:length(tables)) {
	MEAN<-"MEAN"
	SD<-"SD"
	DPtr<-"DPtr"
	TABLEname<-file_path_sans_ext(tables[TABLE])
	print(paste0("Reading table for set ",TABLEname))
	DP<-read.table(tables[TABLE], header=T)
	DP<-DP[complete.cases(DP),]
	summary(DP$DP)
	mean(DP$DP)
	sd(DP$DP)
	png(paste0(TABLEname,"-boxplot.png"))
	boxplot(DP$DP,data=DP, main=TABLEname, cex.main=1, ylab="DP")
	dev.off()
	### Calculate our treshold as follows
	mean<-mean(DP$DP)
	sd<-sd(DP$DP)
	dptr<-mean+5*sd
	MEAN<-rbind(MEAN, mean)
	SD<-rbind(SD, sd)
	DPtr<-rbind(DPtr, dptr)
	write.table(cbind(MEAN,SD,DPtr), paste(TABLEname,".calc",sep=""), sep="\t", col.names=F, row.names=F, quote=F)
}
