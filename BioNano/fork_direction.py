#!/usr/bin/env python
import sys
import numpy as np

def overlap(a, b):
	'''
	Return boolean for whether interval b
	overlaps interval a
	'''
	if (a[0] <= b[0] <= a[1]) or (b[0] <= a[0] <= b[1]):
		return True
	else:
		return False

def parse_dline(line):
	'Format and parse line of d file'
	mysplit = line.split()
	chrom = mysplit[0]
	if chrom == 'chr23':
		chrom = 'chrX'
	if chrom == 'chr24':
		chrom = 'chrY'
	start = int(mysplit[1])
	end = int(mysplit[2])
	sign = mysplit[4]
	strength = float(mysplit[5])
	return chrom, start, end, sign, strength

def fork_score(lf, mlf, rf, mrf):
	'''
	(left forks /mean left forks) - (right forks / mean right forks)
	'''
	score = (lf / mlf) - (rf / mrf)
	return score


def write_fork_scores(dfile, hfile, mfile, ofile, ds):
	'''
	Calculate fork score per bin:
	ds: boolean for direction strength
		if True use sum of direction strength instead 
		of raw fork counts
	'''
	OUT = open(ofile, 'w')
	# Get mean left forks and mean right forks
	with open(mfile, 'r') as m:
		m.readline()
		msplit = m.readline().split()
		mean_lf = float(msplit[0])
		mean_rf = float(msplit[1])

	with open(hfile, 'r') as H: 
		for i, hline in enumerate(H):
			if i%100 == 0:
				print 'On line: ' + str(i)
			leftcnt = 0
			rightcnt = 0
			ds_left_sum = 0
			ds_right_sum = 0
			h_chrom, h_start, h_end = hline.split()
			intervalH = [int(h_start), int(h_end)]
			with open(dfile, 'r') as D:
				for j, dline in enumerate(D):
					d_chrom, d_start, d_end, d_sign, d_strength = parse_dline(dline)
					intervalD = [d_start, d_end]
					# Check if segment overlaps current bin
					if (h_chrom == d_chrom) and overlap(intervalD, intervalH):
						# Add fork
						if d_sign == '+':
							if ds:
								ds_right_sum += d_strength
							else:
								rightcnt += 1
						elif d_sign == '-':
							if ds:
								ds_left_sum += d_strength
							else:
								leftcnt +=1
						else:
							print 'ERROR: something is fucky'
							sys.exit()
			if ds:
				fs = fork_score(ds_left_sum, mean_lf, ds_right_sum, mean_rf)
			else:
				fs = fork_score(leftcnt, mean_lf, rightcnt, mean_rf)
			# Write bin and fork_score to .bedGraph
			OUT.write('\t'.join([h_chrom, h_start, h_end, str(fs)]) + '\n')

	OUT.close()


def main():

	dfile = 'direction_multi_segment_merge_clean_sorted.bed'
	hfile = 'hg19_100kb_sorted.bed'
	# mfile = 'fork_mean.txt'
	mfile = 'fork_ds_means.txt'
	# ofile = 'fork_direction.bedGraph'
	ofile = 'fork_direction_ds.bedGraph'

	# write_fork_scores(dfile, hfile, mfile, ofile)
	ds = True
	write_fork_scores(dfile, hfile, mfile, ofile, ds)






	
if __name__ == '__main__':
	main()

