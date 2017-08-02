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
	return chrom, start, end, sign


def get_mean_forks(dfile, hfile):
	'''
	Calculate mean number of left forks and right 
	forks per 100kb bin
	'''
	leftforks = []
	rightforks = []

	with open(hfile, 'r') as H: 
		for i, hline in enumerate(H):
			if i%100 == 0:
				print 'On line: ' + str(i)
			leftcnt = 0
			rightcnt = 0
			h_chrom, h_start, h_end = hline.split()
			intervalH = [int(h_start), int(h_end)]
			with open(dfile, 'r') as D:
				for j, dline in enumerate(D):
					d_chrom, d_start, d_end, d_sign = parse_dline(dline)
					intervalD = [d_start, d_end]
					# Check if segment overlaps current bin
					if (h_chrom == d_chrom) and overlap(intervalD, intervalH):
						# Add fork
						if d_sign == '+':
							rightcnt += 1
						elif d_sign == '-':
							leftcnt +=1
						else:
							print 'ERROR: something is fucky'
							sys.exit()
			leftforks.append(leftcnt)
			rightforks.append(rightcnt)

	leftmean = np.mean(np.array(leftforks))
	rightmean = np.mean(np.array(rightforks))
	print 'Length of leftforks array: ' + str(len(leftforks)) 
	print 'Length of rightforks array: ' + str(len(rightforks))
	return leftmean, rightmean

def main():

	dfile = 'direction_multi_segment_merge_clean_sorted.bed'
	hfile = 'hg19_100kb_sorted.bed'

	l_mean, r_mean = get_mean_forks(dfile, hfile)
	with open('fork_means.txt', 'w') as t:
		t.write('mean_leftforks\tmean_rightforks\n')
		t.write(str(l_mean) + '\t' + str(r_mean) + '\n')

if __name__ == '__main__':
	main()

