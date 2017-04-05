#!/usr/bin/env python

import argparse
import re
import numpy as np
parser = argparse.ArgumentParser(description='create bedfile for bionano data')
args = parser.parse_args()


def main():

	# Open bnx file with red signal data
	counter = 1
	OUT = open('red_labels.bed', 'w')
	BNX = open('Molecules (647dU).bnx', 'r')
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
			if len(red_line.split()) > 2:

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

					for lab_idx in range(len(aligned_red)):
						OUT.write(chrom + '\t' + str(aligned_red[lab_idx]) + 
						'\t' + str(aligned_red[lab_idx] + 1) + '\t' +  intensity[lab_idx] + '\n')

					# Check if xmap entries ever have reversed refpos
					if RefStartPos >= float(xline[6]):
						print 'WTF'
						quit()
				XMAP.close()				
		counter += 1


if __name__ == '__main__':
	main()