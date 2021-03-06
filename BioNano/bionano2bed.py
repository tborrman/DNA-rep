#!/usr/bin/env python

import argparse
import re
import numpy as np
parser = argparse.ArgumentParser(description='create bedfile for bionano data')
parser.add_argument('-b', help='input bnx file', type=str, required=True)
parser.add_argument('-o', help='output bed file', type=str, required=True)
parser.add_argument('-e', help='exclude singleton labels', type=bool, default=False)
args = parser.parse_args()


def filter_red_labels(red_labels):
	'''
	Get index of red labels to keep
	'''
	idx_to_keep = []
	for idx in range(len(red_labels)):
		total_neighbors = 0
		test_label = red_labels[idx]
		for query_label in red_labels:
			if abs(test_label - query_label) < 20000 :
				total_neighbors +=1
		if total_neighbors > 5:
			idx_to_keep.append(idx)
	return idx_to_keep


def main():

	# Open bnx file with red signal data
	counter = 1
	OUT = open(args.o, 'w')
	BNX = open(args.b, 'r')
	if args.e:
		min_redlabel = 3
	else:
		min_redlabel = 2
	for line in BNX:
		if counter % 1000 == 0:
			print 'On line ' + str(counter)
		if line[0] == '0':
			# Save molecule ID
			mol_id = line.split()[1]
			red_line = next(BNX)
			if red_line[0] != '1':
				print 'ERROR in line format'
				quit()
			if len(red_line.split()) > min_redlabel:

				XMAP = open('all.xmap', 'r')
				xmap_str = XMAP.read()
				srch_obj =  re.search('\n\d+\t'+mol_id+'\t.+', xmap_str)
				if srch_obj:
					xline = srch_obj.group().split()
					chrom = 'chr' + xline[2]
					RefStartPos = float(xline[5])
					RefEndPos = float(xline[6])
					Orientation = xline[7]

					# Get green label info
					green_line = next(BNX)
					if green_line[0] != '2':
						print 'ERROR in line format'
						quit()
					if Orientation == '+':
						first_green = float(green_line.split()[1])
						alignment_factor = RefStartPos - first_green
						red_array = np.array(map(float,red_line.split()))
						aligned_red = red_array[1:-1] + alignment_factor
					
					elif Orientation == '-':
						first_green = float(green_line.split()[1])
						red_array = np.array(map(float,red_line.split()))[1:-1]
						align_steps = red_array - first_green
						aligned_red_bwd= RefEndPos - align_steps
						aligned_red = aligned_red_bwd[::-1]
					else:
						'ERROR no orientation'
						quit()

					idx_red_to_keep = filter_red_labels(aligned_red)
					filtered_aligned_red = [aligned_red[i] for i in idx_red_to_keep]

					# Get label intensity
					next(BNX)
					next(BNX)
					intensity_line = next(BNX)
					if intensity_line.split()[0] != 'QX12':
						print 'ERROR reading intensity'
						quit()
					intensity = intensity_line.split()[1:]

					if Orientation == '-':
						intensity.reverse()

					filtered_intensity = [intensity[i] for i in idx_red_to_keep]

					for lab_idx in range(len(filtered_aligned_red)):
						OUT.write(chrom + '\t' + str(filtered_aligned_red[lab_idx]) + 
						'\t' + str(filtered_aligned_red[lab_idx] + 1) + '\t' +  filtered_intensity[lab_idx] + '\n')

					# Check if xmap entries ever have reversed refpos
					if RefStartPos >= float(xline[6]):
						print 'WTF'
						quit()
				XMAP.close()				
		counter += 1


if __name__ == '__main__':
	main()