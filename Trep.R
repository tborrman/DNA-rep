#!/usr/bin/env Rscript

#################################################
# A script to plot Treps vs genomic position
#
# REQUIRES:
# ori_data_2.0.txt in working directory
#
# ARGUMENTS:
# timecourse_data (ex. chrIV_all_timepts.txt):
#	- tab-delimited text file
#	- col 1: chromosome
# 	- col 2: bin start
# 	- col 3: bin end
# 	- col 4-14: copy number at time point
#
# USAGE:
# ./Trep.R [timecourse_data] 
#
# EXAMPLE:
# ./Trep.R chrIV_all_timepts.txt
#
##################################################

################## BEGIN PROGRAM #################

# Get command line ARGUMENTS
args <- commandArgs()
timecourse_file <- args[6]

# Sigmoid function with adjustable parameters
sigmoid <- function(params, minutes) {
	a = params[1]
	b = params[2]
	c = params[3]
	d = params[4]
	return(a + (d/(1 + exp(-b*(minutes -c)))))
}

# Get time course data
wkDir <- getwd()
time_course <- read.table(paste(wkDir, timecourse_file, sep="/"), header=FALSE, sep="\t", comment.char="")
minutes <- c(5,10,15,20,25,30,35,40,45,60,90)

# Get origin data
raw_origins <- read.table(paste(wkDir, "ori_data_2.0.txt", sep="/"), header=TRUE, sep="\t", comment.char="")
# Select only CONFIRMED origins
confirmed_origins <- raw_origins[which(raw_origins$status == "Confirmed"),]
# Select for chromosome
chr <- substr(as.character(time_course[1,1]), 4, nchar(as.character(time_course[1,1])))
chr <- as.numeric(as.roman(chr))
chr_confirmed_origins <- confirmed_origins[which(confirmed_origins$chr == chr),]

# Get Treps for all 1 kb bins of chromosome
Treps <- c()
genomic_pos <- c()

for (i in 1:nrow(time_course)) {
	if (i%%100 == 0) {
		print(paste("Estimating parameters for row: ", i, sep=""))
	}
	genomic_pos <- c(genomic_pos, time_course[i, 2] + 500)
	copy_number <- as.numeric(time_course[i,4:ncol(time_course)])
	# Get non-linear least square estimates for parameters of sigmoid function
	fitModel <- nls(copy_number ~  a+(d/(1 + exp(-b*(minutes-c)))), start=list(a=1,b=.5, c=40, d=1), trace=FALSE, control = nls.control(warnOnly = TRUE))
	params <- coef(fitModel)
	# Save Treps
	Trep <- params[3]
	# Set upper bound of 60 minutes for Treps to handle estimates that didn't converge (can change later)
	if (Trep > 60) {
		Trep <- 60
	}
	Treps <- c(Treps, Trep)
}

# Get origin positions
origin_pos <- chr_confirmed_origins$center
origin_pos_kb <- origin_pos/1000
# Plot Treps vs genomic position
genomic_pos_kb <- genomic_pos/1000
plotFile <- paste("Treps_vs_genomic_pos_chr_", chr, ".pdf", sep="")
plotPath <- paste(wkDir, plotFile, sep="/")
pdf(plotPath, width=20, height=5)
	par(mar=c(5,5,4,2) + 0.1)
	title <- paste("Chromosome ", as.roman(chr), sep="")
	plot(genomic_pos_kb, Treps, main=title, pch=20, xlab="Genomic position (kb)", ylab= "Trep (minutes)", ylim=c(60,10), xlim=c(0, genomic_pos_kb[length(genomic_pos_kb)]),cex.lab=1.5,cex.main=1.5)
	axis(1, at=seq(0, genomic_pos_kb[length(genomic_pos_kb)], by = 100))
	par(xpd=TRUE)
	points(origin_pos_kb, rep(62,length(origin_pos_kb)), pch=21, bg="white", cex=2, lwd=2)
dev.off();




