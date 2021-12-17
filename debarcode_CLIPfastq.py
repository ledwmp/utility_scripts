import sys
"""Arguments: debarcody.py fastq_file barcode_file
Splits a fastq depending on barcodes defined in barcode_file
"""

class barcode_writer:
	"""Object to store a barcode and open file that fastq entry can write to.
	"""
	def __init__(self,base_fh: str, barcode: str, barcode_number: str):
		self._basefh = base_fh.split("/")[-1]
		self._barcode = barcode
		self._barcode_number = barcode_number
		self._barcode_counter = 0
		self.open_file()
	def open_file(self):
		"""Opens file numbered with barcode_number.
		"""
		self._file_handle = self._barcode_number + "_" + self._basefh
		self._open_file = open(self._file_handle,"a")
	def write_fastq(self,fastq: tuple):
		for i in fastq:
			self._open_file.write(i+"\n")
		self._barcode_counter += 1
	def close_file(self):
		"""Closes barcode number.
		Returns:
			Number of barcodes: int
		"""
		self._open_file.close()
		return self._barcode_counter

class fastq_parser:
	"""Takes CL arguments <fastq_file>.fastq and <barcode_file>.txt with structure:
	barcode_RT_orientatin\tbarcode_number\tbarcode_RC
	and creates barcode_writer objects. Then iterates through fastq feeding
	matching fastqs to barcode_writer objects.
	"""
	def __init__(self,fastq_fh: str, barcode_fh: str,pair: bool):
		self._fastq_fh = fastq_fh
		self._barcode_fh = barcode_fh
		self._pair = pair
		self.populate_barcodes()
		self.iterate_fastq()
		self.close_barcodes()
	def populate_barcodes(self):
		"""Opens <barcode_file.txt> and populates one dictionary with
		barcode:barcode_number pairs and second dictionary with
		barcode:barcode_object pairs.
		"""
		self._sample_index = {}
		with open(self) as r:
			for line in r:
				new_line = line.split("\t")
				self._sample_index[new_line[2].strip()] = new_line[1].strip()
		r.close()
		self._sample_index["bad_barcode"] = "bad_barcode"
		self._sample_barcodes = {j:barcode_writer(i,j) \
		 							for i,j in self._sample_index.items()}
	def iterate_fastq(self):
		"""Walks through fastq and passes entries to respective barcode objects.
		"""
		#is interleaved paired-end fastq
		if self._pair:
			with open(self._fastq_fh) as r:
				for line in r:
					read_name = line.strip()
					read = next(r).strip()
					qual_name = next(r).strip()
					qual = next(r).strip()
					current_bc = read[5:9]
					current_umi = read[0:5]+read[9:11]
					read_name_2 = next(r).strip()
					read_2 = next(r).strip()
					qual_name_2 = next(r).strip()
					qual_2 = next(r).strip()
					if current_bc in self._sample_index.keys():
						self._sample_barcodes[self._sample_index[current_bc]].write_fastq((\
							read_name.split(" ")[0]+":umi:"+current_umi+" "+read_name.split(" ")[1],\
							read[11:],\
							qual_name,\
							qual[11:],\
							read_name_2.split(" ")[0]+":umi:"+current_umi+" "+read_name_2.split(" ")[1],\
							read_2,\
							qual_name_2,\
							qual_2,\
								))
					else:
						self._sample_barcodes["bad_barcode"].write_fastq((\
							read_name.split(" ")[0]+":umi:"+current_umi+" "+read_name.split(" ")[1],\
							read[11:],\
							qual_name,\
							qual[11:],\
							read_name_2.split(" ")[0]+":umi:"+current_umi+" "+read_name_2.split(" ")[1],\
							read_2,\
							qual_name_2,\
							qual_2,\
								))
		#is single end fastq
		else:
			with open(self._fastq_fh) as r:
				for line in r:
					read_name = line.strip()
					read = next(r).strip()
					qual_name = next(r).strip()
					qual = next(r).strip()
					current_bc = read[5:9]
					current_umi = read[0:5]+read[9:11]
					if current_bc in self._sample_index.keys():
						self._sample_barcodes[self._sample_index[current_bc]].write_fastq((\
							read_name.split(" ")[0]+":umi:"+current_umi+" "+read_name.split(" ")[1],\
							read[11:],\
							qual_name,\
							qual[11:],\
								))
					else:
						self._sample_barcodes["bad_barcode"].write_fastq((\
							read_name.split(" ")[0]+":umi:"+current_umi+" "+read_name.split(" ")[1],\
							read[11:],\
							qual_name,\
							qual[11:],\
								))
	def close_barcodes(self):
		"""Calls close_file method on barcode objects and dumps log file.
		"""
		tmp_dict = {j._barcode_number:j.close_file() \
										for i,j in self._sample_barcodes.items()}
		self._logfile = self._fastq_fh.split(".fastq")[0] + "_debarcode.log"
		with open(self._logfile,"a") as r:
			for i,j in tmp_dict.items():
				r.write(i+"\t"+j+"\n")
def main():
	barcode_file = sys.argv[-1]
	fastq_files = sys.argv[1:-1]
	for i in fastq_files:
		parser = fastq_parser(i,barcode_file,True)

if __name__ == "__main__":
	main()
