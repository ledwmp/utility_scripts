#!/usr/bin python
import sys
import numpy as np
from scipy.stats import chi2_contingency,fisher_exact
"""
arguments input_norm_nf.py
		CLIP.count.txt
		input.count.txt
		CLIP.flagstat.txt
		input.flagstat.txt
        CLIP.input_norm.bed
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

def parse_fc_lines(l_clip: str, l_input: str, CLIP_reads: int, input_reads: int):
	"""Parses lines from featureCounts and outputs a tuple that can be output to
	final bed file:
	Args:
		l_clip (str): raw fc line
		l_input (str): raw fc line
		CLIP_reads (int): number of CLIP reads in dataset
		input_reads (int): number of input reads in dataset
	Returns:
		(chrom, \
		start, \
		stop, \
		namepeak, \
		readCLIP, \
		readinput, \
		pvalue, \
		test statistic, \
		Fisher/Chi2, \
		neglog10pvalue, \
		log2FC, \
		entropy, \
		strand, \
		)
	"""
	#ENSG00000284662	1	685679	686673	-	995	0
	l_clip_split = l_clip.strip().split("\t")
	l_input_split = l_input.strip().split("\t")
	if l_clip_split[0:3] != l_input_split[0:3]:
		raise RuntimeError( \
					"featureCount files do not match between CLIP and input" \
					)
	chrom = l_clip_split[1]
	namepeak = l_clip_split[0]
	start = str(int(l_clip_split[2])-1)
	#start = l_clip_split[2]
	stop = l_clip_split[3]
	strand = l_clip_split[4]

	a = int(l_input_split[-1])
	c = int(l_clip_split[-1])
	b = input_reads - a
	d = CLIP_reads - c
	"""
	2by2 contingency table:
						    peak notpeak
					input	|a	|b
					CLIP	|c	|d
	"""
	#normalize reads by cpm
	norm_CLIP = c / CLIP_reads * 1.e6
	norm_input = a / input_reads * 1.e6
	if a >= 1:
		#fold-change
		fold_over_input = norm_CLIP/norm_input
		#from biorxiv paper, IPi * log2(IPi/Qi) where
		#IPi = fraction IP reads in region and Qi = fraction input reads in region
		entropy = (float(c)/CLIP_reads)* \
						np.log2((c/CLIP_reads)/(a/input_reads))
	else:
		#fc,throws a 1 in denominator if no input reads detected in peak
		norm_input = (a+1) / input_reads * 1.e6
		fold_over_input = norm_CLIP/norm_input
		entropy = (float(c)/CLIP_reads)* \
						np.log2((c/CLIP_reads)/((a+1)/input_reads))
	log2fc = np.log2(fold_over_input)
	contin_array = [[a,b],[c,d]]
	if a>=1 and c>=1:
		#chi2 test with Yates' correction
		statistic,pvalue,dof,exp = chi2_contingency( \
										contin_array, correction=True)
		test_type = "CHI2"
	else:
		#fisher's exact test, two-sided
		statistic,pvalue = fisher_exact(contin_array)
		test_type = "FISHER"
	if pvalue == 0.0:
		neglog10pvalue = 400.0
	else:
		neglog10pvalue = -1.0*np.log10(pvalue)
	return (chrom, \
			start, \
			stop, \
			namepeak, \
			str(c), \
			str(a), \
			str(pvalue), \
			str(statistic), \
			test_type, \
			str(neglog10pvalue), \
			str(log2fc), \
			str(entropy), \
			strand, \
			)

def main():
	CLIP_reads = fetch_flagstat_reads(sys.argv[3])
	input_reads = fetch_flagstat_reads(sys.argv[4])
	out_file = open(sys.argv[5],"a")
	with open(sys.argv[1]) as CLIP_file, open(sys.argv[2]) as input_file:
		#skip header in featureCounts files
		next(CLIP_file)
		next(input_file)
		next(CLIP_file)
		next(input_file)
		for line_clip,line_input in zip(CLIP_file,input_file):
			out = parse_fc_lines(line_clip,line_input,CLIP_reads,input_reads)
			out_file.write("\t".join(out)+"\n")
	out_file.close()

if __name__ == "__main__":
	main()
