# A script to find correlations in origin data

library(ggplot2);

# Set paths
begPath <- "/Users/User/Research/DNArep";
raw_ori <-read.table(paste(begPath, "/Data/ori_data_1.3.txt", sep=""), header=TRUE, sep="\t", comment.char="");
# Remove columns we will ignore in correlation 
remove <- c(4:47,52, 59:62, 64, 65, 67:ncol(raw_ori));
prep_ori <- raw_ori[-remove];
# Sum the macapline count
macalpine_MCM_count = prep_ori$rhind_belsky_MacAlpine.Upstream.MCM.Count + prep_ori$rhind_belsky_MacAlpine.Downstream.MCM.Count;
prep_ori <- cbind(prep_ori, macalpine_MCM_count);
belsky_rhind_MCM_count = prep_ori$rhind_belsky_Rhind.Upstream.MCM.Count + prep_ori$rhind_belsky_Rhind.Downstream.MCM.Count;
prep_ori <- cbind(prep_ori, belsky_rhind_MCM_count);
prep_ori <- prep_ori[-c(4,5,6,7)];
# Rename columns for easier reading
colnames(prep_ori)[6] <- "dsSPD4_MCM_count";
colnames(prep_ori)[9] <- "dsSPD8_MCM_count";

# dsSPD4_MCM_count vs t1/2
dsSPD4 <- prep_ori[c(6,11,10)];
ggplot(dsSPD4, aes(x = nieduszynski_T1_2, y = dsSPD4_MCM_count, color = donaldson_WT_yku70_diff_plus3))+
         geom_point(size=3, shape = 16) + #represent the data with points
         scale_colour_gradient(low="blue3", high="yellow") +
         labs(x = "Nieduszynski T 1/2", y="dsSPD4 MCM count", colour = "(WT - yku70) + 3")
         
# dsSPD8_MCM_count vs t1/2
dsSPD8 <- prep_ori[c(9,11,10)];
ggplot(dsSPD8, aes(x = nieduszynski_T1_2, y = dsSPD8_MCM_count, color = donaldson_WT_yku70_diff_plus3))+
  geom_point(size=3, shape = 16) + #represent the data with points
  scale_colour_gradient(low="blue3", high="yellow") +
  labs(x = "Nieduszynski T 1/2", y="dsSPD8 MCM count", colour = "(WT - yku70) + 3")

# macalpine_MCM_count vs t1/2
macalpine <- prep_ori[c(12,11,10)];
ggplot(macalpine, aes(x = nieduszynski_T1_2, y = macalpine_MCM_count, color = donaldson_WT_yku70_diff_plus3))+
  geom_point(size=3, shape = 16) + #represent the data with points
  scale_colour_gradient(low="blue3", high="yellow") +
  labs(x = "Nieduszynski T 1/2", y="Macalpine MCM count", colour = "(WT - yku70) + 3")      

# belsky_rhind_MCM_count vs t1/2
belsky <- prep_ori[c(13,11,10)];
ggplot(belsky, aes(x = nieduszynski_T1_2, y = belsky_rhind_MCM_count, color = donaldson_WT_yku70_diff_plus3))+
  geom_point(size=3, shape = 16) + #represent the data with points
  scale_colour_gradient(low="blue3", high="yellow") +
  labs(x = "Nieduszynski T 1/2", y="Belsky's Rhind MCM count", colour = "(WT - yku70) + 3")       
         
       
# Compute correlation matrix
prep_ori <- prep_ori[4:ncol(prep_ori)];
cor_matrix <- cor(prep_ori, use="pairwise.complete.obs", method="pearson");
# Plot correlations
#pairs(prep_ori);
# Write results
write.table(cor_matrix, paste(begPath, "/Results/MCM/MCM_T1_2_correlations.txt", sep=""), sep='\t', row.names=TRUE, col.names=NA);