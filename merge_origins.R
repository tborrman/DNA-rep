# A script to merge replication origin data from multiple datasets.

# Set paths
begPath <- "/Users/User/Research/DNArep";
wkDir <- paste(begPath, "/Data", sep="");

# Functions
source(paste(begPath, '/Code/functions_data_processing.R', sep=""));

# READ IN DATA ------------------------------------------------------------

# Read in skeletal replication origin data from OriDB
skeleton <- read.table(paste(wkDir, "/oriDB_origins_S_cer.txt", sep=""), header=TRUE, sep="\t", comment.char="", quote="");
# Read in yabuki_raghu.txt to get time of replication data
Trep <- read.table(paste(wkDir, "/yabuki_raghu.txt", sep=""), header=TRUE, sep="\t", comment.char="", quote="");
# Read in Whitehouse efficiency data
whitehouse <- read.table(paste(wkDir, "/whitehouse/whitehouse.txt", sep=""), header=TRUE, sep="\t", comment.char="", quote="");
# Read in Chen efficiency data
chen <- read.table(paste(wkDir, "/chen.txt", sep=""), header=TRUE, sep="\t", comment.char="", quote="");


# VERSION 1.0 -------------------------------------------------------------



# Aggregate data
# Initialize variables
dirty_data <- c();

# Give unique ID to each origin
ID <- seq(1:nrow(skeleton));
skeleton <- cbind(ID, skeleton);

# Rename columns to clarify datasets that were aggregated
colnames(Trep) <- paste("Trep", colnames(Trep), sep="_");
colnames(whitehouse) <- paste("whitehouse", colnames(whitehouse), sep="_");
colnames(chen) <- paste("chen", colnames(chen), sep="_");

# Add OriCenter to OriDB data
center <- ceiling(((skeleton$end - skeleton$start)/2) + skeleton$start);
skeleton <- cbind(skeleton[1:2], center, skeleton[3:ncol(skeleton)]);

# Format chromosome numbers in Chen data
chen_chr <- c();
for(i in 1:nrow(chen)){
  if(chen$chen_ChrName[i] == "chrI") {
    chen_chr[i] <- 1;
  } else if (chen$chen_ChrName[i] == "chrII") {
    chen_chr[i] <- 2;
  } else if (chen$chen_ChrName[i] == "chrIII") {
    chen_chr[i] <- 3;
  } else if (chen$chen_ChrName[i] == "chrIV") {
    chen_chr[i] <- 4;
  } else if (chen$chen_ChrName[i] == "chrV") {
    chen_chr[i] <- 5;
  } else if (chen$chen_ChrName[i] == "chrVI") {
    chen_chr[i] <- 6;
  } else if (chen$chen_ChrName[i] == "chrVII") {
    chen_chr[i] <- 7;
  } else if (chen$chen_ChrName[i] == "chrVIII") {
    chen_chr[i] <- 8;
  } else if (chen$chen_ChrName[i] == "chrIX") {
    chen_chr[i] <- 9;
  } else if (chen$chen_ChrName[i] == "chrX") {
    chen_chr[i] <- 10;
  } else if (chen$chen_ChrName[i] == "chrXI") {
    chen_chr[i] <- 11;
  } else if (chen$chen_ChrName[i] == "chrXII") {
    chen_chr[i] <- 12;
  } else if (chen$chen_ChrName[i] == "chrXIII") {
    chen_chr[i] <- 13;
  } else if (chen$chen_ChrName[i] == "chrXIV") {
    chen_chr[i] <- 14;
  } else if (chen$chen_ChrName[i] == "chrXV") {
    chen_chr[i] <- 15;
  } else if (chen$chen_ChrName[i] == "chrXVI") {
    chen_chr[i] <- 16;
  }
}
chen <- cbind(chen_chr, chen);

# Since not all origins have names, mapping by bp position requires breaking up data by chromosomes.
totalChrom <- 16;
for (chrom in c(1:totalChrom)) {
  skeleton_block <- skeleton[which(skeleton$chr == chrom),];
  

# STEP 1 ------------------------------------------------------------------


  # Merge yabuki raghu trep data by starting position
  Trep_block <- Trep[which(Trep$Trep_chr == chrom),];
  # Need to match to closest bp from OriDB
  nearest_bps <- sapply(Trep_block$Trep_start, function(x) which.min(abs(x-skeleton_block$start)));
  # Add matching bp column for merging by bp
  start <- skeleton_block$start[nearest_bps];
  Trep_match <- cbind(start, Trep_block);
  # Check matching start locations for duplicates and remove duplicates with greatest origin difference
  duplicates <- which(duplicated(Trep_match$start));
  if(length(duplicates) != 0) {
    duplicates <- c((duplicates[1]-1), duplicates);
    keep <- which.min(abs(Trep_match[duplicates, "start"] - Trep_match[duplicates, "Trep_start"]));
    remove <- duplicates[-keep];
    Trep_noDups <- Trep_match[-remove, ];
  } else { 
    Trep_noDups <- Trep_match;
  }
  # Remove any data with bp difference greater than 6000
  mismatches <- which(abs(Trep_noDups$start - Trep_noDups$Trep_start) > 6000);
  
  if (length(mismatches) != 0) {
    Trep_clean <- Trep_noDups[-mismatches,];
  } else {
    Trep_clean <- Trep_noDups;
  }
  data_1_block <- merge.with.order(skeleton_block, Trep_clean, by="start", all.x = TRUE, sort=F ,keep_order = 1);
  
  # Merge whitehouse data by starting position
  whitehouse_block <- whitehouse[which(whitehouse$whitehouse_chr == chrom),];
  # Need to match to closest bp from OriDB
  nearest_bps <- sapply(whitehouse_block$whitehouse_WT_A_bp, function(x) which.min(abs(x-skeleton_block$center)));
  # Add matching bp column for merging by bp
  center <- skeleton_block$center[nearest_bps];
  whitehouse_match <- cbind(center, whitehouse_block);
  # Check matching start locations for duplicates and remove duplicates with greatest origin difference
  duplicates <- which(duplicated(whitehouse_match$center));
  if(length(duplicates) != 0) {
    duplicates <- c((duplicates[1]-1), duplicates);
    keep <- which.min(abs(whitehouse_match[duplicates, "center"] - whitehouse_match[duplicates, "whitehouse_WT_A_bp"]));
    remove <- duplicates[-keep];
    whitehouse_noDups <- whitehouse_match[-remove, ];
   } else { 
     whitehouse_noDups <- whitehouse_match;
   }
  # Remove any data with bp difference greater than 6000
  mismatches <- which(abs(whitehouse_noDups$center - whitehouse_noDups$whitehouse_WT_A_bp) > 6000);
  if (length(mismatches) != 0) {
    whitehouse_clean <- whitehouse_noDups[-mismatches,];
  } else {
    whitehouse_clean <- whitehouse_noDups;
  }
  # Merge data
  data_2_block <- merge.with.order(data_1_block, whitehouse_clean, by = "center", all.x=TRUE, sort=F, keep_order = 1);
  
  # Merge chen data by center 
  chen_block <- chen[which(chen$chen_chr == chrom),];
  # Need to match to closest bp from OriDB
  nearest_bps <- sapply(chen_block$chen_OriCent, function(x) which.min(abs(x-skeleton_block$center)));
  # Add matching bp column for merging by bp
  center <- skeleton_block$center[nearest_bps];
  chen_match <- cbind(center, chen_block);
  # Check matching start locations for duplicates and remove duplicates with greatest origin difference
  duplicates <- which(duplicated(chen_match$center));
  if(length(duplicates) != 0) {
    duplicates <- c((duplicates[1]-1), duplicates);
    keep <- which.min(abs(chen_match[duplicates, "center"] - chen_match[duplicates, "chen_OriCent"]));
    remove <- duplicates[-keep];
    chen_noDups <- chen_match[-remove, ];
  } else { 
    chen_noDups <- chen_match;
  }
  # Remove any data with bp difference greater than 6000
  mismatches <- which(abs(chen_noDups$center - chen_noDups$chen_OriCent) > 6000);
  if (length(mismatches) != 0) {
    chen_clean <- chen_noDups[-mismatches,];
  } else {
    chen_clean <- chen_noDups;
  }
  # Merge data
  data_3_block <- merge.with.order(data_2_block, chen_clean, by = "center", all.x=TRUE, sort=F, keep_order = 1);
  
  # Bind the block data for each chromosome
  dirty_data <- rbind(dirty_data, data_3_block);
}

# STEP 2 ------------------------------------------------------------------


# Clean the data
# Rearrange order of for ID, chr, start, center, end columns
clean_data <- dirty_data
colnames(clean_data) [1] <- "ID";
clean_data[1] <- dirty_data$ID;
colnames(clean_data)[2] <- "chr";
clean_data[2] <- dirty_data$chr;
colnames(clean_data)[3] <- "start";
clean_data[3] <- dirty_data$start;
colnames(clean_data)[4] <- "center";
clean_data[4] <- dirty_data$center;
colnames(clean_data)[5] <- "end";
clean_data[5] <- dirty_data$end

# Check matching names (a similar check was performed for chen and Trep resulting in zero name discrepancies)
# We show there exists name discrepancies in whitehouse data which will be manually adjusted for
for (i in 1:nrow(clean_data)) {
  if(is.na(clean_data$whitehouse_Name[i])) {
    #print(clean_data$whitehouse_Name[i]);
    
  }
  else if (as.character(clean_data$whitehouse_Name[i]) == as.character(clean_data$name[i])) {
    #print(as.character(clean_data$whitehouse_Name[i]));
  } 
  else {
    #print(as.character(clean_data$whitehouse_Name[i]));
  }
}

# Adjust whitehouse data for names: ARS434, ARS1022, ARS1218, ARS1628
clean_data[189, 20:28] <- clean_data[188, 20:28];
clean_data[188, 20:28] <- NA;
clean_data[454, 20:28] <- clean_data[453, 20:28];
clean_data[453, 20:28] <- NA;
clean_data[538, 20:28] <- clean_data[537, 20:28];
clean_data[537, 20:28] <- NA;
clean_data[820, 20:28] <- clean_data[821, 20:28];
clean_data[821, 20:28] <- NA;


# Output the aggregate table
write.table(clean_data, paste(wkDir, "/ori_data_1.0.txt", sep=""), sep='\t', row.names=F);


# Run bp check to see if centers of origins greatly differ between datasets.
# Check bp difference between Ori and whitehouse_A experiment
bp_diff <- abs(clean_data$center - clean_data$whitehouse_WT_A_bp);
mean_whitehouseA_diff <- mean(bp_diff, na.rm=TRUE);
median_whitehouseA_diff <- median(bp_diff, na.rm=TRUE);
std_whitehouseA_diff <- sd(bp_diff, na.rm=TRUE);
max_whitehouseA_diff <- max(bp_diff, na.rm=TRUE);

bp_difference <- na.omit(bp_diff);

hist(bp_difference, breaks = 80);

bp_distance <- abs(clean_data$whitehouse_WT_A_bp - clean_data$center);
OEM <- clean_data$whitehouse_WT_A_efficiency;

plot(bp_distance, OEM);







