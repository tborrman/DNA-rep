#!/usr/bin/env python
import argparse
import numpy as np

parser = argparse.ArgumentParser(description= 'Calculate fork direction profile from the replication data')
parser.add_argument('-b', help = 'binned genome bed (ex. hg19_100kb_sorted.bed)', type=str, required=True)
parser.add_argument('-r', help = 'replication timing bedGraph (ex. RT_HeLaS3_hg19_sorted.bedGraph)', type=str, required=True)
parser.add_argument('-o', help = 'output file (ex. rep2forks.bedGraph)', type=str, required=True)
args = parser.parse_args()


def overlap(a, b):
	'''
	Return boolean for whether interval b
	overlaps interval a
	'''
	if (a[0] <= b[0] <= a[1]) or (b[0] <= a[0] <= b[1]):
		return True
	else:
		return False


def parse_rline(line):
	'Format and parse line of d file'
	mysplit = line.split()
	chrom = mysplit[0]
	if chrom == 'chr23':
		chrom = 'chrX'
	if chrom == 'chr24':
		chrom = 'chrY'
	start = int(mysplit[1])
	end = int(mysplit[2])
	RT = float(mysplit[3])
	return chrom, start, end, RT

def fork_score(up, down):
	mean_dwn = np.mean(np.array(down))
	mean_up = np.mean(np.array(up))
	fs = mean_dwn - mean_up
	return fs


def rep2fork(hfile, rfile, ofile):
	'''
	Calculate a fork direction profile
	from the replication timing data
	'''
	window = 200000 # Look window size upstream and downstream of bin
	OUT = open(ofile, 'w')
	with open(hfile, 'r') as h:
		for i, hline in enumerate(h):
			if i%100 == 0:
				print 'On line: ' + str(i)
			h_chrom, h_start, h_end = hline.split()
			intervalH = [int(h_start) - window, int(h_end) + window]
			upstream = []
			downstream = []
			with open(rfile, 'r') as r:
				for j, rline in enumerate(r):
					r_chrom, r_start, r_end, r_RT = parse_rline(rline)
					intervalR = [r_start, r_end]
					# Check if RT region overlaps current bin
					if (h_chrom == r_chrom) and overlap(intervalR, intervalH):
						if r_end <= int(h_start) :
							upstream.append(r_RT)
						elif r_start >= int(h_end):
							downstream.append(r_RT)
			if len(upstream) == 0 or len(downstream) == 0:
				OUT.write('\t'.join([h_chrom, h_start, h_end, 'NA']) + '\n')
			else:
				fs = fork_score(upstream, downstream)
				OUT.write('\t'.join([h_chrom, h_start, h_end, str(fs)]) + '\n')


def main():

	rep2fork(args.b, args.r, args.o)



if __name__ == '__main__':
	main()


