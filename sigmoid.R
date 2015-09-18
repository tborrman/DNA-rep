#!/usr/bin/env Rscript

################################################
# A script to make sigmoid plot for selected bin
# in DNA replication time course data
# 
# ARGUMENTS:
# datafile (ex. chrIV_all_timepts.txt):
#	- tab-delimited text file
#	- col 1: chromosome
# 	- col 2: bin start
# 	- col 3: bin end
# 	- col 4-14: copy number at time point
# start (ex. 1000):
#	- bp start of bin
#	- col 2 of datafile
# 
# USAGE:
# ./sigmoid.R [datafile] [start]
# 
# EXAMPLE:
# ./sigmoid.R chrIV_all_timepts.txt 5000
#
##############################################	

############## BEGIN PROGRAM #################

# Get command line ARGUMENTS
args <- commandArgs()
datafile <- args[6]
start <- as.numeric(args[7])

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
data <- read.table(paste(wkDir, datafile, sep="/"), header=FALSE, sep="\t", comment.char="")
minutes <- c(5,10,15,20,25,30,35,40,45,60,90)

# Get copy number data for user-defined bin
colnames(data)[2] <- "bin_start"
bin_data <- data[which(data$bin_start == start),]
copy_number <- as.numeric(bin_data[,4:ncol(bin_data)])
chr <- bin_data[1,1]
if (nrow(bin_data) == 0) {
	print("ERROR: bin start missing from table")
	quit()
}
print(copy_number)
# Get non-linear least square estimates for parameters of sigmoid function
fitModel <- nls(copy_number ~  a+(d/(1 + exp(-b*(minutes-c)))), start=list(a=1,b=.5, c=40, d=1), trace=TRUE, control = nls.control(warnOnly = TRUE))
params <- coef(fitModel)
print(params)
# Make list of minute values for function plotting
new_minutes <- seq(0,95);
# Get new model fit values
fit_values <- sigmoid(params, new_minutes);
# Plot result
profileFile <- paste("sigmoidfit_", start, ".pdf", sep="");
profilePath <- paste(wkDir, profileFile, sep="/");
end <- start + 1000
pdf(profilePath, width=6, height=6);
    title <- paste("Copy number (blue) and fit (red)\n for ", chr, " and bin ", start, "-", end , sep="");
    plot(minutes, copy_number, main = title,  ylab= "Copy number", xlab="Minutes", pch = 20, col = "blue");
    lines(new_minutes, fit_values, pch=20, col="red");
dev.off();

# Print Trep
Trep <- params[3]
print(paste("Trep: ", Trep, sep=""))
