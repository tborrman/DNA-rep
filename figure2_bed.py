#!/usr/bin/env python




def idx_start_and_end_segments(label_pos, idx_to_keep):

	segment_idx = []
	begin_segment = True

	for i in range(len(idx_to_keep)-1):

		current_pos = label_pos[idx_to_keep[i]]
		next_pos = label_pos[idx_to_keep[i+1]]

		if abs(current_pos - next_pos) < 30000 and begin_segment:
			segment = [idx_to_keep[i]]
			begin_segment=False
		elif abs(current_pos - next_pos) >= 30000 and not begin_segment:
			segment.append(idx_to_keep[i])
			begin_segment = True
			segment_idx.append(segment)
		# Last entry
		if i + 1 == len(idx_to_keep) -1 and not begin_segment:
			segment.append(idx_to_keep[i+1])
			segment_idx.append(segment)


	return segment_idx



def filter_red_labels(red_labels):
	'''
	Get index of red labels to keep
	'''
	idx_to_keep = []
	for idx in range(len(red_labels)):
		total_neighbors = 0
		test_label = red_labels[idx]
		for query_label in red_labels:
			if abs(test_label - query_label) < 50000 :
				total_neighbors +=1
		if total_neighbors > 4:
			idx_to_keep.append(idx)
	return idx_to_keep


def main():
	for chrom in map(str,range(1 ,24)):

		IN = open('H9_sync2_ch1_labels_chr' + chrom+'.txt', 'r')
		OUT = open('H9_sync2_ch1_labels_chr' + chrom+'.bed', 'w')
		myhash = {}
		for line in IN:
			splitline = line.split()
			if splitline[0] in myhash:
				myhash[splitline[0]].append(float(splitline[1]))
			else:
				myhash[splitline[0]] = [float(splitline[1])]

		for mol_id in myhash:
			positions = myhash[mol_id]
			
			idx_red_to_keep= filter_red_labels(positions)
			if len(idx_red_to_keep) > 4:
				segment_idx = idx_start_and_end_segments(positions, idx_red_to_keep)

				for segment_start_end in segment_idx:
					# Get start and end of segment
					segment_start_bp = positions[segment_start_end[0]]
					segment_end_bp = positions[segment_start_end[1]]

					# Write to bedfile
					OUT.write('chr' + chrom +'\t'+str(segment_start_bp)+'\t'+str(segment_end_bp)+'\n')
		IN.close()
		OUT.close()



if __name__ == '__main__':
	main()


