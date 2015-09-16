# Histogram for fragment lengths (dsSPD8a, and dsSPD8b)

begPath <- "/Users/User/Research/DNArep";
fragment_lengths_dsSPD8a <- scan(paste(begPath, "/Data/rhind/MCM/dsSPD8a_mapped.sam_fragment_lengths.txt", sep=""),);
hist(fragment_lengths_dsSPD8a, breaks=250);
fragment_lengths_dsSPD8b <- scan(paste(begPath, "/Data/rhind/MCM/dsSPD8b_mapped.sam_fragment_lengths.txt", sep=""),);
hist(fragment_lengths_dsSPD8b, breaks=250);
