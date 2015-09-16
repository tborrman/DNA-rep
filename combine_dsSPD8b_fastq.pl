# A script to combine fastq files for dsSPD8b pair end reads in rhind data

# Initialize ouput file for single read data
open(OUT1, ">/home/borrmant/scratch/rhind/dsSPD8b_1.fastq");

# Cycle through fastq files
for $i(1..17) {
    if($i < 10) {
	open(IN1, "/home/borrmant/scratch/rhind/Fastq_files/DS_SPD8_NoIndex_L006_R1_00$i.fastq");
	while($_ = <IN1>) {
	    print OUT1 $_;
	} 
    } else {
	open(IN1, "/home/borrmant/scratch/rhind/Fastq_files/DS_SPD8_NoIndex_L006_R1_0$i.fastq");
	while($_ = <IN1>) {
	    print OUT1 $_;
	}
    }
    close IN1;
}
close OUT1;

# Initialize ouput file for mate read data
open(OUT2, ">/home/borrmant/scratch/rhind/dsSPD8b_2.fastq");

# Cycle through fastq files
for $i(1..17) {
    if($i < 10) {
        open(IN2, "/home/borrmant/scratch/rhind/Fastq_files/DS_SPD8_NoIndex_L006_R2_00$i.fastq");
        while($_ = <IN2>) {
            print OUT2 $_;
        }
    } else {
        open(IN2, "/home/borrmant/scratch/rhind/Fastq_files/DS_SPD8_NoIndex_L006_R2_0$i.fastq");
        while($_ = <IN2>) {
            print OUT2 $_;
        }
    }
    close IN2;
}
close OUT2;

    
		  
		  

