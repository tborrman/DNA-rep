# A script to merge ori_data_1.5 with Knott's RPD3 data outputting ori_data_1.6

# Set paths
begPath <- "/Users/User/Research/DNArep";
wkDir <- paste(begPath, "/Data", sep="");

# Functions
source(paste(begPath, '/Code/functions_data_processing.R', sep=""));

# READ IN DATA ------------------------------------------------------------
# Read in ori_data_1.5
ori_data_1.5 <- read.table(paste(wkDir, "/ori_data_1.5.txt", sep=""), header=TRUE, sep="\t", comment.char="");
# Read in RPD3 data
knott <- read.table(paste(wkDir, "/knott/Table_S3.txt", sep=""), header=TRUE, sep="\t", comment.char="");

# Intialize Ouput
ori_data_1.6 <- c();

# Rename columns to identify dataset that is aggregated
colnames(knott) <- paste("knott", colnames(knott), sep="_");

# Create extra column for binary classification of timing difference (1 = significant difference, 0 = no difference)
knott_sig_diff <- rep(1, nrow(knott));
knott <- cbind(knott, knott_sig_diff);

# MERGE ORIGINS
# We will merge by start of ori_data and knott
# Since not all origins have names, mapping by bp position requires breaking up data by chromosomes.
totalChrom <- 16;
for (chrom in c(1:totalChrom)) {
  ori_data_1.5_block <- ori_data_1.5[which(ori_data_1.5$chr == chrom),];
  knott_block <- knott[which(knott$knott_Chrome == chrom),];
  # Merge knott data by bp location
  # Need to match to closest bp from OriDB
  #yang_pos <- yang_block$yang_ori_pos * 1000;
  nearest_bps <- sapply(knott_block$knott_Start, function(x) which.min(abs(x-ori_data_1.5_block$start)));
  # Add matching bp column for merging by bp
  start <- ori_data_1.5_block$start[nearest_bps];
  knott_match <- cbind(start, knott_block);
  
  # Check matching center locations for duplicates and throw error if so
  duplicates <- which(duplicated(knott_match$start));
  if(length(duplicates) != 0) {
    print("ERROR: duplicate origin matches for acs motif start");
    #print("resolving by selecting only closest match of duplicates...");
    #mate_duplicates <- duplicates -1;
    #all_duplicates <- union(mate_duplicates, duplicates);
    #
    #name_col <- 7;
    #ars_names <- as.character(MCM_counts_match[all_duplicates, name_col]);
    # Match back to ori data
    #new_centers <- c();
    #for (name in ars_names) {
     #   new_center <- which(ori_data_1.0_block$name == name);
      #  new_centers <- c(new_centers, new_center);
    #}
    # Replace incorrect centers with new centers
    #ori_center_col <- 4;
    #mcm_center_col <- 1;
    #print(as.character(MCM_counts_match[all_duplicates, 7]));
    #MCM_counts_match[all_duplicates, mcm_center_col] <- ori_data_1.0_block[new_centers, ori_center_col];
  }
  
  # Perform the merge by center and acs_motif
  ori_data_1.6_block <- merge.with.order(ori_data_1.5_block, knott_match, by="start", all.x = TRUE, sort=F ,keep_order = 1);
  
  # Bind merged data for each chromosome together
  ori_data_1.6 <- rbind(ori_data_1.6, ori_data_1.6_block);

}


# Check to see if all ars names match in new dataset
#mismatch_names <- which(as.character(ori_data_1.6$name) == as.character(ori_data_1.6$knott_ARS.Name));

                        
ori_name <- as.character(ori_data_1.6$name);
knott_name <- as.character(ori_data_1.6$knott_ARS.Name);
knott_name[which(is.na(knott_name))] <- "";

ori_othernames <- as.character(ori_data_1.6$othernames);
knott_othernames <- as.character(ori_data_1.6$knott_Alternative.Name);
#knott_othernames[which(is.na(knott_othernames))] <- "";

# Look at all mismatches 
mismatch_names <- which(ori_name != knott_name);
print(paste("There are a total of ", length(mismatch_names)," mismatched names", sep=""));

# Find mismatches where mcm_name is given but ori_name is missing
ori_missing <- intersect(which(ori_name != ""), which(knott_name == ""));
print(paste("There are a total of ", length(ori_missing), " mismatched names where the knott_name alone is missing", sep=""));

mismatch_othernames <- which(ori_othernames != knott_othernames);
print(length(mismatch_othernames));

ori_other_missing <- intersect(which(ori_othernames != ""), which(knott_othernames == ""));
print(length(ori_other_missing));

# After above checking and visual inspection of discrepances merge was successful. 
# Put start column back in appropriate place
ori_data_1.6_done <- ori_data_1.6
colnames(ori_data_1.6_done) [1] <- "ID";
ori_data_1.6_done[1] <- ori_data_1.6$ID;
colnames(ori_data_1.6_done)[2] <- "chr";
ori_data_1.6_done[2] <- ori_data_1.6$chr;
colnames(ori_data_1.6_done)[3] <- "start";
ori_data_1.6_done[3] <- ori_data_1.6$start;
colnames(ori_data_1.6_done)[4] <- "center";
ori_data_1.6_done[4] <- ori_data_1.6$center;
ori_data_1.6 <- ori_data_1.6_done;

# Make all NA entries of sig_difference = 0
ori_data_1.6[which(is.na(ori_data_1.6$knott_sig_diff)), "knott_sig_diff"] <- 0;


# Output the aggregate table
write.table(ori_data_1.6, paste(wkDir, "/ori_data_1.6.txt", sep=""), sep='\t', row.names=F);







