#!/usr/bin/env python
import subprocess

for idx in range(129):
	subprocess.call("bsub -q long -W 24:00 -R 'rusage[mem=5000]' -o " + str(idx) + ".out -e " +str(idx) + ".err " +
		"./direction_bed_seg_new_format.py -b Molecules_split_bnx_" + "{:03d}".format(idx) + " -o Molecules_split_bnx_" + "{:03d}".format(idx) + 
		".bed", shell= True)
