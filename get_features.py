#!/usr/bin/env python



def pixel_width(molID, num_pixles):
	molLength = None
	midBNX = open('93620B2ecli_2017-01-03_15_25_middleSide.bnx', 'r')
	for bnxline in midBNX:
		splitline = bnxline.split()
		if splitline[0] == '0' and splitline[1] == molID:
			molLength = float(splitline[2])
			break
	midBNX.close()
	if not molLength:
		'ERROR: can not find molecule in BNX file'
		quit()
	else:
		pixWidth = molLength / num_pixles
		return pixWidth

def segment_indices(seg):
	idx_list = []
	segFlag = False
	old_level = seg[0]

	for i, level in enumerate(seg[1:-1]):
		idx = i + 1
		if not segFlag and level > old_level:
			seg_idx = []
			segFlag = True
			seg_idx.append(idx)
			old_level = level
		elif segFlag and level < old_level:
			seg_Flag = False
			seg_idx.append(idx - 1)
			old_level = level
			idx_list.append(seg_idx)
	if len(idx_list) == 0:
		print 'No segment in this molecule'
	return idx_list


def get_nick1_distance(molID):
	IN = open('93620B2ecli_2017-01-03_15_25_greenFiltered.bnx', 'r')
	for line in IN:
		if line[0] == '0':
			if line.split()[1] == molID:
				nextline=next(IN)
				IN.close()
				return float(nextline.split()[1])

def get_ref_nick_start(molID):
	IN = open('93620B2ecli_2017-01-03_15_25_greenFiltered_errF_contig1.xmap', 'r')
	for line in IN:
		splitline = line.split()
		if splitline[1] == molID:
			return float(splitline[5])








		

	


	



def main():

	TRACE = open('93620B2ecli_2017-01-03_15_25_middle_traceFile.txt', 'r')
	OUT = open('new_test.txt', 'w')

	TRACE.readline()
	counter = 0
	fiber_num = 0

	for pixelline in TRACE:
			counter +=1
			if counter % 10 == 0:
				print 'On line: ' + str(counter)
			# Skip intensity line
			next(TRACE)
			# Get segmentation level
			segmentline = next(TRACE)
			molecule_ID = pixelline.split()[0]
			if molecule_ID != segmentline.split()[0]:
				print 'ERROR: molecule IDs do not match'
				quit()
			num_pix = len(segmentline.split()) - 1
			pixWidth = pixel_width(molecule_ID, num_pix)
			seg_array = segmentline.split()[1:]
			seg_idx_list = segment_indices(map(float,seg_array))
			
			# Get distance of first nick to subtract from reference start
			nick1_distance = get_nick1_distance(molecule_ID)
			# Map to reference
			ref_nick_start = get_ref_nick_start(molecule_ID)

			# For preserving info for multiple segments in one fiber
			# if ref_nick_start:
			# 	ref_start = ref_nick_start - nick1_distance
			# 	if ref_start > 0 :
			# 		for i, seg_index in enumerate(seg_idx_list):	
			# 			segment_start = ref_start + (pixWidth * (seg_index[0]+1))
			# 			segment_end = ref_start + (pixWidth * (seg_index[1]+1))
			# 			if segment_end - segment_start >= 10000:
			# 				if i == 0:
			# 					fiber_num +=1
			# 				OUT.write(str(fiber_num) + '\t')
			# 				OUT.write(str(segment_start)+'\t'+str(segment_end)+'\n')
			if ref_nick_start:
				ref_start = ref_nick_start - nick1_distance
				if ref_start > 0 :
					for seg_index in seg_idx_list:	
						segment_start = ref_start + (pixWidth * (seg_index[0]+1))
						segment_end = ref_start + (pixWidth * (seg_index[1]+1))
						if segment_end - segment_start >= 10000:
							OUT.write(str(segment_start)+'\t'+str(segment_end)+'\n')

			if counter == 100:
				TRACE.close()
				OUT.close()
				quit()

	TRACE.close()
	OUT.close()




if __name__ == '__main__':
	main()

