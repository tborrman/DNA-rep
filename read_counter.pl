#!/usr/bin/perl

# A script to count the number of overlapping reads for a given window size around a list of given origins
# Input: 
# 1. file with list of origin locations (belsky.txt)
# 2. window size (even integer in bp)
# 3. prefix for sorted  bam file (dsSPD4a, etc. note that an indexed bam file must be located in same directory)
# 4. prefix for associated bam file (dsSPD4b_1, etc.)
# 5. prefix for associated bam file (dsSPD8a, etc.)
# 6. prefix for associated bam file (dsSPD8b, etc.)
# Example: read_counter.pl belsky.txt 1000 dsSPD4a dsSPD4b_1 dsSPD8a dsSPD8b 

# Prerequisites:
# samtools must be installed and binary in PATH

use strict;
use warnings;


# Arguments
my $origin_locs=$ARGV[0];
my $window_size=$ARGV[1];
my $read_map1 = $ARGV[2];
my $read_map2 = $ARGV[3];
my $read_map3 = $ARGV[4];
my $read_map4 = $ARGV[5];

###BEGIN PROGRAM########################################################################
# Create string to represent summed MCM count over datasets (just use first 6 characters)
$read_map1=~ /.{6}/;
my $read_map_1_2 = $&;
$read_map3=~ /.{6}/;
my $read_map_3_4 = $&;
# Assign mapped read files
my $map_file1 = join("", $read_map1, "_mapped_sorted.bam");
my $map_file2 = join("", $read_map2, "_mapped_sorted.bam");
my $map_file3 = join("", $read_map3, "_mapped_sorted.bam");
my $map_file4 = join("", $read_map4, "_mapped_sorted.bam");
# Open location file and file to print results
open(ORIGINS, $origin_locs);
open(OUT, join("", ">rhind_borrman", "_MCM_counts.txt"));

# Set variables for origin location input file (need to be adjusted depending on file)
my $chr_col = 1;
my $acs_loc_col = 3;

# Skip the header line of location file
my $header =  <ORIGINS>;
$header =~ /.+/;

# Add new read count columns to output file 
my $new_header =  join("\t", $&, "${read_map1}_MCM_count", "${read_map2}_MCM_count", "${read_map_1_2}_MCM_count", "${read_map3}_MCM_count", "${read_map4}_MCM_count", "${read_map_3_4}_MCM_count");
print OUT  "$new_header\n"; 
# Cycle through origin locations outputting read count for each location
while ($_ = <ORIGINS>) {
    my @line = split;
    # Extract origin location info
    my $chr = $line[$chr_col];
    my $acs_loc = $line[$acs_loc_col];
    # Get read counts overlapping window around acs motif start location
    my $window_start = $acs_loc - ($window_size/2);
    my $window_end = $acs_loc + ($window_size/2);
    my $mcm_count1 = `samtools view -c $map_file1 chr$chr:$window_start-$window_end`;
    $mcm_count1 =~ /.+/;
    $mcm_count1 = $&;
    my $mcm_count2 = `samtools view -c $map_file2 chr$chr:$window_start-$window_end`;
    $mcm_count2 =~ /.+/;
    $mcm_count2 = $&;
    my $mcm_count_1_2 = $mcm_count1 + $mcm_count2;
    my $mcm_count3 = `samtools view -c $map_file3 chr$chr:$window_start-$window_end`;
    $mcm_count3 =~ /.+/;
    $mcm_count3 = $&;
    my $mcm_count4 = `samtools view -c $map_file4 chr$chr:$window_start-$window_end`;
    $mcm_count4 =~ /.+/;
    $mcm_count4 = $&;
    my $mcm_count_3_4 = $mcm_count3 + $mcm_count4;
    
    # Join results
    my $mcm_results = join("\t", $mcm_count1, $mcm_count2, $mcm_count_1_2, $mcm_count3, $mcm_count4, $mcm_count_3_4);
    # Print mcm_count to file
    push(@line, $mcm_results);
    my $new_line  = join("\t", @line);
    print OUT "$new_line\n";
}

close OUT;


