import pysam
import sys
import numpy as np


item = sys.argv[1]
bam_in = pysam.AlignmentFile(item,"rb")
dm6_path = item.split(".bam")[0]+".dm6.bam"
hg38_path = item.split(".bam")[0]+".hg38.bam"
bam_dm6 = pysam.AlignmentFile(dm6_path,"wb",template=bam_in)
bam_hg38 = pysam.AlignmentFile(hg38_path,"wb",template=bam_in)
dm6_list = []
hg38_list = []
for read in bam_in:
	alignment_score = read.get_tag('AS')
	sequence = read.query_sequence
	chromosome = bam_in.get_reference_name(read.reference_id)
	start = read.reference_start
	if "dm6_" in chromosome:
		bam_dm6.write(read)
		#dm6_list.append(read.query_name)
	else:
		bam_hg38.write(read)
		#hg38_list.append(read.query_name)

hg38_set = set(hg38_list)
dm6_set = set(dm6_list)
print(len(hg38_set))
print(len(dm6_set))
print(len(hg38_set | dm6_set))
