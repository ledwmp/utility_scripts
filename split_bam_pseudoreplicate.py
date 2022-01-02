import pysam
import sys


item = sys.argv[1]
bam_in = pysam.AlignmentFile(item,"rb")
rep1_path = "_".join(item.split("_")[0:2])+"rep1_"+"_".join(item.split("_")[1:])
rep2_path = "_".join(item.split("_")[0:2])+"rep2_"+"_".join(item.split("_")[1:])
bam_rep1 = pysam.AlignmentFile(rep1_path,"wb",template=bam_in)
bam_rep2 = pysam.AlignmentFile(rep2_path,"wb",template=bam_in)

#A00589:346:HL5KMDRXY:1:1103:14407:6684:umi:ACAATTC
for read in bam_in:
	read_name = read.query_name
      #read numbered based on y-unit
	read_number = int(read_name.split(":")[6])
	if read_number % 2 == 1:
		bam_rep1.write(read)
	else:
		bam_rep2.write(read)
bam_in.close()
bam_rep1.close()
bam_rep2.close()

