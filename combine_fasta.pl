# A script to merge chromosome fasta files for cerevisiae_fasta

# Initialize output file
open(OUT, ">/home/borrmant/nearline/DNArep/Data/cerevisiae_genome.fasta");

# Cycle through chromosome files
for $i(1..16) {
    open(IN, "/home/borrmant/nearline/DNArep/Data/cerevisiae_fasta/chr$i.fa");
    while($_ = <IN>) {
	print OUT $_;
    }
    close IN;
}
close OUT;
