#!/usr/bin/env python
import argparse
import sys
import pandas as pd
import numpy as np

parser=argparse.ArgumentParser(description='build 2x2 contingency table')
parser.parse_args()

def rep2dict(fname):
	'''
	Create dictionary of predicted fork direction based on rep
	'''
	d = {}
	with open(fname, 'r') as f:
		for line in f:
			splitline = line.split()
			loc = splitline[0]+':'+splitline[1] 
			fs = splitline[3]
			if fs != 'NA':
				if float(fs) < 0:
					# polarity not sign
					fs = '+'
				elif float(fs) >= 0:
					fs = '-'
				else:
					print 'WTF'
					sys.exit()
			d[loc] = fs
	return d

def nearest_floor(x):
	return int(np.floor(x/100000.0) * 100000)

def test(row, r):
	TL = 0
	FL = 0
	FR = 0
	TR = 0
	chrom = row[0]
	start = str(nearest_floor(row[1]))
	pred_sign = row[4]
	try:
		true_sign = r[chrom+':'+start]
	except KeyError:
		return [TL, FL, FR, TR]
	if pred_sign == '+':
		if pred_sign == true_sign:
			TR = 1
		elif true_sign == 'NA':
			pass
		else:
			FR = 1
	if pred_sign == '-':
		if pred_sign == true_sign:
			TL = 1
		elif true_sign == 'NA':
			pass
		else:
			FL = 1
	return ([TL, FL, FR, TR])

def contingency_table(df, r):
	'''
	 	   Rep Left   Rep Right
		  +---------+----------+
Bio Left  | True L  | False L  |
	      +---------+----------+
Bio Right | False R | True R   |
		  +---------+----------+  
	'''
	OUT = open('prediction_table.txt', 'w')
	OUT.write('\t'.join(['strength>', 'True_left', 'False_left', 'False_right', 'True_right']) + '\n')

	# Strenth tolerance to filter segments
	tols = np.arange(0, max(df.strength), max(df.strength)/10.0)
	for tol in tols:
	 	df_t = df.loc[df['strength'] >= tol]		
		pred_arr = np.array(list(df_t.apply(test, axis=1, args=(r,))))
		# Remove rows with NAs where predictions couldn't be made
		x = pred_arr[~np.all(pred_arr == 0, axis=1),:]
		results = np.sum(x, axis=0)
		OUT.write(str(tol) + '\t'+'\t'.join(map(str,results)) + '\n')
	return


def main():
	# segment file processed from bionano data
	seg_file = 'direction_multi_segment_merge_clean_sorted.bed'
	# fork orientation file predicted directly from replication timing data
	rep_file = 'rep2forks.bedGraph'
	rep_fork = rep2dict(rep_file)

	sdf  = pd.read_csv(seg_file, sep='\t', header=None, 
		names=['chrom', 'start', 'end', 'molid', 'sign', 'strength', 'sigsum'])

	# create contingency table
	contingency_table(sdf, rep_fork)


if __name__ == '__main__':
	main()
