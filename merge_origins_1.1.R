# A script to merge ori_data_1.0 with rhind_borrman_MCM_counts.txt outputting ori_data_1.1

# Set paths
begPath <- "/Users/User/Research/DNArep";
wkDir <- paste(begPath, "/Data", sep="");

# Functions
source(paste(begPath, '/Code/functions_data_processing.R', sep=""));

# READ IN DATA ------------------------------------------------------------
# Read in ori_data_1.0
ori_data_1.0 <- read.table(paste(wkDir, "/ori_data_1.0.txt", sep=""), header=TRUE, sep="\t", comment.char="");
# Read in MCM count data
MCM_counts <- read.table(paste(wkDir, "/rhind/MCM/rhind_borrman_MCM_counts.txt", sep=""), header=TRUE, sep="\t", comment.char="");

# Intialize Ouput
ori_data_1.1 <- c();

# Rename columns to identify dataset that is aggregated
colnames(MCM_counts) <- paste("rhind_belsky", colnames(MCM_counts), sep="_");

# ARS602 has two identical copies in MCM dataset, delete one of them.
ARS602_copies <- which(MCM_counts$rhind_belsky_ars_name == "ARS602");
MCM_counts <- MCM_counts[-ARS602_copies[1],];

# MERGE ORIGINS
# We will merge by the center of origins for ori_data_1.0 and the acs_motif_start for MCM_counts.
# Since not all origins have names, mapping by bp position requires breaking up data by chromosomes.
totalChrom <- 16;
for (chrom in c(1:totalChrom)) {
  ori_data_1.0_block <- ori_data_1.0[which(ori_data_1.0$chr == chrom),];
  MCM_counts_block <- MCM_counts[which(MCM_counts$rhind_belsky_chr == chrom),];
  # Merge MCM count data by ACS motif start position
  # Need to match to closest bp from OriDB
  nearest_bps <- sapply(MCM_counts_block$rhind_belsky_acs_motif_start, function(x) which.min(abs(x-ori_data_1.0_block$center)));
  # Add matching bp column for merging by bp
  center <- ori_data_1.0_block$center[nearest_bps];
  MCM_counts_match <- cbind(center, MCM_counts_block);
  
  # Check matching center locations for duplicates and throw error if so
  duplicates <- which(duplicated(MCM_counts_match$center));
  if(length(duplicates) != 0) {
    print("ERROR: duplicate origin matches for acs motif start");
    print("resolving by matching by ars name...");
    mate_duplicates <- duplicates -1;
    all_duplicates <- union(mate_duplicates, duplicates);
    # Get ars name from MCM count data
    name_col <- 7;
    ars_names <- as.character(MCM_counts_match[all_duplicates, name_col]);
    # Match back to ori data
    new_centers <- c();
    for (name in ars_names) {
        new_center <- which(ori_data_1.0_block$name == name);
        new_centers <- c(new_centers, new_center);
    }
    # Replace incorrect centers with new centers
    ori_center_col <- 4;
    mcm_center_col <- 1;
    print(as.character(MCM_counts_match[all_duplicates, 7]));
    MCM_counts_match[all_duplicates, mcm_center_col] <- ori_data_1.0_block[new_centers, ori_center_col];
  }
  # Perform the merge by center and acs_motif
  ori_data_1.1_block <- merge.with.order(ori_data_1.0_block, MCM_counts_match, by="center", all.x = TRUE, sort=F ,keep_order = 1);
  
  # Bind merged data for each chromosome together
  ori_data_1.1 <- rbind(ori_data_1.1, ori_data_1.1_block);

}

# Check to see if all ars names match in new dataset
#mismatch_names <- which(as.character(ori_data_1.1$name) == as.character(ori_data_1.1$rhind_belsky_ars_name));
                        
ori_name <- as.character(ori_data_1.1$name);
mcm_name <- as.character(ori_data_1.1$rhind_belsky_ars_name);
mcm_name[which(is.na(mcm_name))] <- "";

# Look at all mismatches 
mismatch_names <- which(ori_name != mcm_name);
print(paste("There are a total of ", length(mismatch_names)," mismatched names", sep=""));

# Find mismatches where mcm_name is given but ori_name is missing
ori_missing <- intersect(which(ori_name != ""), which(mcm_name == ""));
print(paste("There are a total of ", length(ori_missing), " mismatched names where the mcm_name alone is missing", sep=""));

# All mismatches are due to mcm_data having no name where ori_data does which is acceptable. 

# Put center column back in appropriate place
ori_data_1.1_done <- ori_data_1.1
colnames(ori_data_1.1_done) [1] <- "ID";
ori_data_1.1_done[1] <- ori_data_1.1$ID;
colnames(ori_data_1.1_done)[2] <- "chr";
ori_data_1.1_done[2] <- ori_data_1.1$chr;
colnames(ori_data_1.1_done)[3] <- "start";
ori_data_1.1_done[3] <- ori_data_1.1$start;
colnames(ori_data_1.1_done)[4] <- "center";
ori_data_1.1_done[4] <- ori_data_1.1$center;
ori_data_1.1 <- ori_data_1.1_done;

# Output the aggregate table
write.table(ori_data_1.1, paste(wkDir, "/ori_data_1.1.txt", sep=""), sep='\t', row.names=F);







