# A script to find correlations in origin data

# Set paths
begPath <- "/Users/User/Research/DNArep";
raw_ori <-read.table(paste(begPath, "/Data/ori_data_1.1.txt", sep=""), header=TRUE, sep="\t", comment.char="");
# Remove columns we will ignore in correlation 
remove <- c(1:16, 18, 20:21, 23:37, 39:46, 50:54, 56:57);
prep_ori <- raw_ori[-remove];
# Sum the macapline count
macalpine_MCM_count = prep_ori$rhind_belsky_MacAlpine.Upstream.MCM.Count + prep_ori$rhind_belsky_MacAlpine.Downstream.MCM.Count;
prep_ori <- cbind(prep_ori, macalpine_MCM_count);
prep_ori <- prep_ori[-c(6,7)];
# Rename columns for easier reading
colnames(prep_ori)[4] <- "chen_WT_efficiency";
colnames(prep_ori)[5] <- "belsky_efficiency";
colnames(prep_ori)[6] <- "dsSPD4_MCM_count";
colnames(prep_ori)[7] <- "dsSPD8_MCM_count";
# Compute correlation matrix
cor_matrix <- cor(prep_ori, use="pairwise.complete.obs", method="pearson");
# Plot correlations
pairs(prep_ori);
# Write results
write.table(cor_matrix, paste(begPath, "/Results/MCM/efficiency_mcm_corr_matrix.txt", sep=""), sep='\t', row.names=TRUE, col.names=NA);