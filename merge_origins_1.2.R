# A script to merge ori_data_1.1 with KU Data from Donaldson outputting ori_data_1.2

# Set paths
begPath <- "/Users/User/Research/DNArep";
wkDir <- paste(begPath, "/Data", sep="");

# Functions
source(paste(begPath, '/Code/functions_data_processing.R', sep=""));

# Initialize output
ori_data_1.2 <- c();

# READ IN DATA ------------------------------------------------------------
# Read in ori_data_1.1
ori_data_1.1 <- read.table(paste(wkDir, "/ori_data_1.1.txt", sep=""), header=TRUE, sep="\t", comment.char="");
# Read in Ku Data by chromosome
chrom_number <- 1:16;
for (chrom in chrom_number) {
  Ku_data <- read.table(paste(wkDir, "/donaldson/KU_data/", chrom, "meanTrepNR.txt", sep=""), header=FALSE, sep="\t", comment.char="");
  colnames(Ku_data) <- c("donaldson_kb", "B", "C", "D", "E", "donaldson_WT_Trep", "donaldson_yku70_Trep");
  # Remove superfluous columns
  remove <- 2:5;
  Ku_data <- Ku_data[-remove];
  # Add chromosome column
  donaldson_chrom <- rep(chrom, nrow(Ku_data));
  Ku_data <- cbind(donaldson_chrom, Ku_data);
  # Get Trep difference + 3 (yku70 mutant had a 3 minute delay)
  donaldson_WT_yku70_diff_plus3 <- (Ku_data$donaldson_WT_Trep - Ku_data$donaldson_yku70_Trep) + 3;
  Ku_data <- cbind(Ku_data, donaldson_WT_yku70_diff_plus3);
  # Change kb to bp
  Ku_data$donaldson_kb <- Ku_data$donaldson_kb*1000;
  colnames(Ku_data)[2] <- "donaldson_bp";
  # Get chromosome block of ori_data
  ori_data_1.1_block <- ori_data_1.1[which(ori_data_1.1$chr == chrom),];
  # MERGE DATA
  # Wish to find the nearest Trep measurement for each origin by bp location
  nearest_bps <- sapply(ori_data_1.1_block$center, function(x) which.min(abs(x-Ku_data$donaldson_bp)));
  Ku_matched_data <- Ku_data[nearest_bps, ];
  # Merge data by nearest bp
  ori_data_1.2_block <- cbind(ori_data_1.1_block, Ku_matched_data);
  
  # Bind merged data for each chromosome together
  ori_data_1.2 <- rbind(ori_data_1.2, ori_data_1.2_block);
}

# Output the aggregate table
write.table(ori_data_1.2, paste(wkDir, "/ori_data_1.2.txt", sep=""), sep='\t', row.names=F);

