#!/usr/bin python
import sys
import numpy as np
"""
arguments input_norm_nf.py
		CLIP1.count.txt
		input1.count.txt
		CLIP2.count.txt
		input2.count.txt
		CLIP1.flagstat.txt
		input1.flagstat.txt
		CLIP2.flagstat.txt
		input2.flagstat.txt
        CLIP.final.bed
		CLIP.filter_fc.bed
"""
def fetch_flagstat_reads(file_path: str):
	"""Opens flagstat file and returns the number of reads in the bam
	"""
	with open(file_path) as r:
		#48861 + 0 read1
		for line in r:
			if " read1" in line:
				return int(line.split("+")[0].strip())
	raise RuntimeError("File is probably not samtools flagstat output.")

def parse_fc_lines(l_clip1: str, l_input1: str, l_clip2: str, l_input2: str, \
 	   CLIP_reads1: int, input_reads1: int, CLIP_reads2: int, input_reads2: int):
	"""Parses lines from featureCounts and outputs a tuple that can be output to
	final bed file:
	Args:
		l_clip1 (str): raw fc line
		l_input1 (str): raw fc line
		l_clip2 (str): raw fc line
		l_input2 (str): raw fc line
		CLIP1_reads (int): number of CLIP reads in dataset
		input1_reads (int): number of input reads in dataset
		CLIP2_reads (int): number of CLIP reads in dataset
		input2_reads (int): number of input reads in dataset
	Returns:
		(chrom, \
		start, \
		stop, \
		namepeak, \
		log2cpmCLIP1|log2cpmCLIP2|log2cpmCLIPave|log2fc1|log2fc2|log2fcave, \
		strand, \
		)
	"""
	#ENSG00000284662	1	685679	686673	-	995	0
	l_clip_split1 = l_clip1.strip().split("\t")
	l_input_split1 = l_input1.strip().split("\t")
	l_clip_split2 = l_clip2.strip().split("\t")
	l_input_split2 = l_input2.strip().split("\t")
	if not (l_clip_split1[0:3] == l_input_split1[0:3] == \
	 								l_clip_split2[0:3] == l_input_split2[0:3]):
		raise RuntimeError( \
					"featureCount files do not match between CLIP and input" \
					)
	chrom = l_clip_split1[1]
	namepeak = l_clip_split1[0]
	start = str(int(l_clip_split1[2])-1)
	#start = l_clip_split[2]
	stop = l_clip_split1[3]
	strand = l_clip_split1[4]

	peakCLIP1 = int(l_clip_split1[-1])
	peakinput1 = int(l_input_split1[-1])
	peakCLIP2 = int(l_clip_split2[-1])
	peakinput2 = int(l_input_split2[-1])

	#normalize reads by cpm
	norm_CLIP1 = peakCLIP1 / CLIP_reads1 * 1.e6
	norm_input1 = peakinput1 / input_reads1 * 1.e6
	norm_CLIP2 = peakCLIP2 / CLIP_reads2 * 1.e6
	norm_input2 = peakinput2 / input_reads2 * 1.e6

	if norm_input1 >= 1:
		#fold-change
		fold_over_input1 = norm_CLIP1/norm_input1
	else:
		#fc,throws a 1 in denominator if no input reads detected in peak
		norm_input1 = (peakinput1+1) / input_reads1 * 1.e6
		fold_over_input1 = norm_CLIP1/norm_input1

	if norm_input2 >= 1:
		#fold-change
		fold_over_input2 = norm_CLIP2/norm_input2
	else:
		#fc,throws a 1 in denominator if no input reads detected in peak
		norm_input2 = (peakinput2+1) / input_reads2 * 1.e6
		fold_over_input2 = norm_CLIP2/norm_input2

	log2fc1 = str(np.log2(fold_over_input1))
	log2fc2 = str(np.log2(fold_over_input2))
	log2fcave = str(np.log2((fold_over_input1+fold_over_input2)/2.))
	log2cpm1 = str(np.log2(norm_CLIP1))
	log2cpm2 = str(np.log2(norm_CLIP2))
	log2cpmave = str(np.log2((norm_CLIP1+norm_CLIP2)/2.))

	score = "|".join([log2fc1,log2fc2,log2fcave,log2cpm1,log2cpm2,log2cpmave])
	return (chrom, \
			start, \
			stop, \
			namepeak, \
			score, \
			strand, \
			)

def main():
	CLIP1_reads = fetch_flagstat_reads(sys.argv[5])
	input1_reads = fetch_flagstat_reads(sys.argv[6])
	CLIP2_reads = fetch_flagstat_reads(sys.argv[7])
	input2_reads = fetch_flagstat_reads(sys.argv[8])
	out_file_final = open(sys.argv[9],"a")
	out_file_filter = open(sys.argv[10],"a")
	with open(sys.argv[1]) as CLIP1_file, open(sys.argv[2]) as input1_file, \
			open(sys.argv[3]) as CLIP2_file, open(sys.argv[4]) as input2_file:
		#skip header in featureCounts files
		next(CLIP1_file),next(input1_file),next(CLIP2_file),next(input2_file)
		next(CLIP1_file),next(input1_file),next(CLIP2_file),next(input2_file)
		for line_clip1,line_input1,line_clip2,line_input2 \
		 				in zip(CLIP1_file,input1_file,CLIP2_file,input2_file):
			out = parse_fc_lines( \
					line_clip1, \
					line_input1, \
					line_clip2, \
					line_input2, \
					CLIP1_reads, \
					input1_reads, \
					CLIP2_reads, \
					input2_reads,\
					)
			out_file_final.write("\t".join(out)+"\n")
			if float(out[4].split("|")[2]) > 0.32:
				out_file_filter.write("\t".join(out)+"\n")
	out_file_final.close()
	out_file_filter.close()
if __name__ == "__main__":
	main()
