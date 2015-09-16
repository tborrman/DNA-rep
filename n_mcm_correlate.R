# A script to find correlations in origin data

library(ggplot2);

# Set paths
begPath <- "/Users/User/Research/DNArep";
raw_ori <-read.table(paste(begPath, "/Data/ori_data_1.6.txt", sep=""), header=TRUE, sep="\t", comment.char="");

# Remove origins missing parameter n info
ori_data_n <- raw_ori[-which(is.na(raw_ori$yang_n)),];

# Remove rDNA origin skewing MCM ChIP seq
ori_data_clean <- ori_data_n[-which(rownames(ori_data_n) == 534),];

ori_data <- ori_data_clean[c("dsSPD4_MCM_count", "yang_n", "donaldson_WT_yku70_diff_plus3", "knott_sig_diff")];


# For plotting change values of knott_sig_diff
ori_data[which(ori_data$knott_sig_diff == 1), "knott_sig_diff"] <- "Rpd3 dependent";
ori_data[which(ori_data$knott_sig_diff == 0), "knott_sig_diff"] <- "Non-Rpd3 dependent";

# dsSPD4_MCM_count vs n
ggplot(ori_data, aes(x = yang_n, y = dsSPD4_MCM_count, color = donaldson_WT_yku70_diff_plus3))+
         geom_point(size=3, shape = 16) + #represent the data with points
         scale_colour_gradient(low="blue3", high="yellow") +
         labs(x = "Yang parameter n", y="dsSPD4 MCM count", colour = "(WT - yku70) + 3")

ggplot(ori_data, aes(x = yang_n, y = dsSPD4_MCM_count, color = knott_sig_diff))+
  geom_point(size=3, shape = 16) + #represent the data with points
  #scale_colour_gradient(low="blue3", high="yellow") +
  labs(x = "Yang parameter n", y="dsSPD4 MCM count", colour = "Timing Dependence")

# Compute correlation matrix
macalpine_mcm_count <- ori_data_clean$rhind_belsky_MacAlpine.Upstream.MCM.Count + ori_data_clean$rhind_belsky_MacAlpine.Downstream.MCM.Count;
ori_data_clean <- cbind(ori_data_clean, macalpine_mcm_count);
cor_ori <- ori_data_clean[c("yang_n", "dsSPD4a_MCM_count", "dsSPD4b_1_MCM_count", "dsSPD4_MCM_count", "dsSPD8a_MCM_count", "dsSPD8b_MCM_count", "dsSPD8_MCM_count", "macalpine_mcm_count")]
cor_matrix <- cor(cor_ori, use="pairwise.complete.obs", method="pearson");
# Plot correlations
#pairs(prep_ori);
# Write results
write.table(cor_matrix, paste(begPath, "/Results/MCM/MCM_n_correlations.txt", sep=""), sep='\t', row.names=TRUE, col.names=NA);

# Compute correlation matrix
ori_data_filter <- ori_data_clean[-which(ori_data_clean$donaldson_WT_yku70_diff_plus3 > 3),];

cor_ori_filter <- ori_data_filter[c("yang_n", "dsSPD4a_MCM_count", "dsSPD4b_1_MCM_count", "dsSPD4_MCM_count", "dsSPD8a_MCM_count", "dsSPD8b_MCM_count", "dsSPD8_MCM_count", "macalpine_mcm_count")]
cor_matrix_filter <- cor(cor_ori_filter, use="pairwise.complete.obs", method="pearson");
# Plot correlations
#pairs(prep_ori);
# Write results
write.table(cor_matrix_filter, paste(begPath, "/Results/MCM/MCM_n_filtered_correlations.txt", sep=""), sep='\t', row.names=TRUE, col.names=NA);