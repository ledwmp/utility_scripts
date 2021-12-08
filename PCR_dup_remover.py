import sys
import re
"""
arguments PCR_dup_remover.py logfile,file1,file2,file3,...
Removes duplicate umis from SAM file if they have the same read start position
"""

class pcr_collapser:
    """Class to parse sam file and output PCR-deduplicated reads based on umi and
    read start positions.
    """
	def __init__(self, logfile: str, fastq_fh: str, pair: bool):
        self._logfile = logfile
        self._fastq_fh = fastq_fh
        self._pair = pair
        self.open_outfile()
    @staticmethod
    def sam_to_bit(flag: int):
        """Method to convert bitwise sam flags to binary
        """
        return str(bin(flag))[2:].zfill(17)
    @staticmethod
    def parse_cigar(cigar: str):
        """Method to find the number of reference bases consumed by the cigar string,
        to be used for finding read start positions aligning to reverse strand.
        """
        cigars = re.findall('(\d+)([MDN=X])', cigar)
        return sum([int(i[0]) for i in cigars])
    def open_outfile(self):
        """Method that iterates through sam file and writes reads to out file
        if they don't exist in set. Also writes stats to log file.
        """
        self._umis = set()
        self._deduplicated_umis = 0
        self._removed_umis = 0
        self._dedup_fh = self._fastq_fh.split(".fastq") + "_dedup.fastq"
        self._out = open(self._dedup_fh,"a")
        self._in = open(self._fastq_fh,"r")
        for line in self._in:
            if line[0] == "@":
                self._out.write(line)
            else:
                new_line = line.strip().split("\t")
                bit_string = sam_to_bit(int(new_line[1]))
                if bit_string[-3] == 1:
                    continue #unmapped read
                elif self._pair:
                    if bit_string[-2] != 1:
                        continue #if paired and not mapped in proper pair
                    elif bit_string[-8] == 1:
                        continue #if paired and second read
                #mapped forward strand
                if bit_string[-5] != 1:
                    seq_start = new_line[3]
                #mapped reverse strand
                else:
                    seq_start = str(int(new_line[3]) + parse_cigar(new_line[5]))
                umi = new_line[0].split("umi:")[1]
                #umi,chrom,seq_start
                test_fragment = str(umi+"_"+new_line[2]+"_"+seq_start)
                #haven't seen umi
                if test_fragment not in self._umis:
                    self._out.write(line)
                    self._umis.add(test_fragment)
                    self._deduplicated_umis += 1
                #seen umi
                else:
                    self._removed_umis += 1
        self._in.close()
        self._out.close()
        with open(self._logfile,"a") as r:
            r.write(self._fastq_fh.split("/")[-1]+"\n"+\
                "Deduplicated_umis:"+str(self._deduplicated_umis)+"\n"+\
                "Removed_umis:"+str(self._removed_umis)+"\n")
        print("Deduplicated_umis:"+str(deduplicated_umis))
        print("Removed_umis:"+str(removed_umis))

def main():
    for i in sys.argv[2:]:
        dedup = pcr_collapser(sys.argv[1],i,True)

if __name__ == "__main__":
    main()
