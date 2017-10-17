#!/usr/bin/env python
import sys


IN = open(sys.argv[1], 'r')
for line in IN:
	splitline = line.split()
	if float(splitline[1]) > 0 and float(splitline[2]) > 0 :
		start = int(round(float(splitline[1]))) -1
		end = int(round(float(splitline[2])))
		print splitline[0]+'\t'+str(start)+'\t'+ str(end) +'\txname\t0\t'+splitline[4]



