# A script to find correlations in origin data

# Set paths
begPath <- "/cygwin64/home/Tyler/Research/DNArep";
# Read in data 
raw_ori <-read.table(paste(begPath, "/Data/ori_data_1.0.txt", sep=""), header=TRUE, sep="\t", comment.char="");
# Remove columns we will ignore in correlation 
remove <- c(1:15, 20:21, 23, 25:27, 29:35);
prep_ori <- raw_ori[-remove];
# Rename columns for easier reading
colnames(prep_ori) <- c("Trep.Yabuki.min", "Trep.Yabuki.RI", "Trep.Raghu.min", "Trep.Raghu.RI", "whitehouse.eff.A", "whitehouse.eff.B", "whitehouse.Raghu.Trep","chen.whitehouse.A", "chen.whitehouse.B", "chen.WT");
# Compute correlation matrix
cor_matrix <- cor(prep_ori, use="pairwise.complete.obs", method="pearson");
# Plot correlations
pairs(prep_ori);
# Write results
write.table(cor_matrix, paste(begPath, "/Results/ori_1.0_corr_matrix.txt", sep=""), sep='\t', row.names=TRUE, col.names=NA);