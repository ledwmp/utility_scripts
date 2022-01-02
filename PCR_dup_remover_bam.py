import sys
import re
import pysam
"""
arguments PCR_dup_remover.py logfile,file1,file2,file3,...
Removes duplicate umis from SAM file if they have the same read start position
"""

class pcr_collapser:
    """Parses sam file and output PCR-deduplicated reads based on umi and
    read start positions.
    """
    def __init__(self, logfile: str, sam_fh: str, pair: bool):
        self._logfile = logfile
        self._sam_fh = sam_fh
        self._pair = pair
        self.open_outfile()
    @staticmethod
    def sam_to_bit(flag: int):
        """Converts bitwise sam flags to binary
        """
        return str(bin(flag))[2:].zfill(17)
    @staticmethod
    def parse_cigar(cigar: str):
        """Finds the number of reference bases consumed by the cigar string,
        to be used for finding read start positions aligning to reverse strand.
        """
        cigars = re.findall('(\d+)([MDN=X])', cigar)
        return sum([int(i[0]) for i in cigars])
    def open_outfile(self):
        """Iterates through sam file and writes reads to out file
        if they don't exist in set. Also writes stats to log file.
        """
        self._umis = set()
        self._deduplicated_umis = 0
        self._removed_umis = 0
        self._dedup_fh = self._sam_fh.split(".bam")[0] + ".dedup.bam"
        self._bam_in = pysam.AlignmentFile(self._sam_fh,"rb")
        self._bam_out = pysam.AlignmentFile(self._dedup_fh,"wb",template=self._bam_in)
        for read in self._bam_in:
            bit_string = self.sam_to_bit(read.flag)
            if int(bit_string[-3]) == 1:
                continue #unmapped read
            elif self._pair:
                if int(bit_string[-2]) != 1:
                    continue #if paired and not mapped in proper pair
                elif int(bit_string[-8]) == 1:
                    continue #if paired and second read
            #mapped forward strand
            if int(bit_string[-5]) != 1:
                seq_start = read.reference_start
            #mapped reverse strand
            else:
                seq_start = str(int(read.reference_start) + self.parse_cigar(read.cigarstring))
            umi = read.query_name.split("umi:")[1]
            #umi,chrom,seq_start
            test_fragment = str(umi+"_"+str(read.reference_id)+"_"+str(seq_start))
            #haven't seen umi
            if test_fragment not in self._umis:
                self._bam_out.write(read)
                self._umis.add(test_fragment)
                self._deduplicated_umis += 1
            #seen umi
            else:
                self._removed_umis += 1
        self._bam_in.close()
        self._bam_out.close()
        with open(self._logfile,"a") as r:
            r.write(self._sam_fh+"\n"+\
                "Deduplicated_umis:"+str(self._deduplicated_umis)+"\n"+\
                "Removed_umis:"+str(self._removed_umis)+"\n")
        print("Deduplicated_umis:"+str(self._deduplicated_umis))
        print("Removed_umis:"+str(self._removed_umis))

def main():
    for i in sys.argv[2:]:
        dedup = pcr_collapser(sys.argv[1],i,True)

if __name__ == "__main__":
    main()
