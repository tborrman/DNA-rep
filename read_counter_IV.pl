#!/usr/bin/perl

# A script to count the number of overlapping reads for a given window size around a list of given origins
# This script count reads given by single bam file that has been sorted and indexed
# Input: 
# 1. file with list of origin locations (ori_data_1.8.txt)
# 2. window size (even integer in bp)
# 3. prefix for sorted  bam file (labib_c, note that an indexed bam file must be located in same directory)

# Example: read_counter_IV.pl ori_data_1.8.txt 1000 labib_c 

# Prerequisites:
# samtools must be installed and binary in PATH

use strict;
use warnings;


# Arguments
my $origin_locs=$ARGV[0];
my $window_size=$ARGV[1];
my $read_map = $ARGV[2];

# Variables
my $window_start;
my $window_end;

###BEGIN PROGRAM########################################################################
# Assign mapped read files
my $map_file = join("", $read_map, "_sorted.bam");

# Open location file and file to print results
open(ORIGINS, $origin_locs);
open(OUT,  ">ori_data_1.9.txt");

# Set variables for origin location input file (need to be adjusted depending on file)
my $chr_col = 1;
my $acs_loc_col = 3;

# Skip the header line of location file
my $header =  <ORIGINS>;
$header =~ /.+/;

# Add new read count columns to output file 
my $new_header =  join("\t", $&, "${read_map}_MCM_counts");
print OUT  "$new_header\n"; 
# Cycle through origin locations outputting read count for each location
while ($_ = <ORIGINS>) {
    my @line = split;
    # Extract origin location info
    my $chr = $line[$chr_col];
    my $acs_loc = $line[$acs_loc_col];
    # Get read counts overlapping window around acs motif start location
    if($acs_loc < ($window_size/2)) {
	$window_start = 0;
	$window_end = $window_size;
    } else {
	$window_start = $acs_loc - ($window_size/2);
	$window_end = $acs_loc + ($window_size/2);
    }
    my $mcm_count = `samtools view -c $map_file chr$chr:$window_start-$window_end`;
    $mcm_count =~ /.+/;
    $mcm_count = $&;
    
    # Print mcm_count to file
    push(@line, $mcm_count);
    my $new_line  = join("\t", @line);
    print OUT "$new_line\n";
}

close OUT;


