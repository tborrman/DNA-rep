#!/usr/bin/perl
# A script to output vector of fragment lengths for pair end reads
# Input: SAM file
# Example: fragment_length.pl dsSPD8a_mapped.sam
use strict;
use warnings;


my $sam_file = $ARGV[0];
open(SAM, "$sam_file");
# The fragment column in sam files is the 9th column, (use 8 due to 0-index); 
my $chr_col = 2;
my $frag_col = 8;

open(OUT, ">${sam_file}_fragment_lengths.txt");
while($_ = <SAM>) {
    my @read = split;
    
    # Only read lines that were mapped to chromosomes
    if($read[$chr_col]=~/chr[0-9]+/) {
	# Record fragment length
	print OUT "$read[$frag_col]\n";
	# Skip next line which will be the mate for the pair
	$_ = <SAM>;
    }

}

close SAM;
close OUT;    
    
