# Script to output correlation plots between belsky and rhind mcm results using 1kb window
begPath <- "/Users/User/Research/DNArep";
raw_data <-read.table(paste(begPath, "/Data/rhind/MCM/rhind_borrman_MCM_counts.txt", sep=""), header=TRUE, sep="\t", comment.char="");
macalpine_mcm_count <- raw_data[,10] + raw_data[,11];
belsky_rhind_mcm_count <- raw_data[,12] + raw_data[,13];
raw_data <- cbind(raw_data, macalpine_mcm_count);
raw_data <- cbind(raw_data, belsky_rhind_mcm_count);
remove <- seq(1,14);
data <- raw_data[, -remove];

# Compute correlation matrix
cor_matrix <- cor(data, method="pearson");
# Plot correlations
pairs(data);
# Write results
write.table(cor_matrix, paste(begPath, "/Results/mcm_count_correlations.txt", sep=""), sep='\t', row.names=TRUE, col.names=NA);