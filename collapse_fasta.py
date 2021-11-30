from bloom_filter import BloomFilter
import sys

kmer_len = 40

class bloom_object:
	"""
	Object to hold a bloom filter and associated list of fasta names
	"""
	def __init__(self,kmers_contained: tuple):
		self.name_list = []
		self.name_list.append(kmers_contained[0])
		self.kmer_len = len(kmers_contained[1][0])
		self.max_elements = len(kmers_contained[1])*4*self.kmer_len
		self.bloom = BloomFilter(max_elements=self.max_elements,error_rate=1./self.kmer_len)
		for i in kmers_contained[1]:
			self.bloom.add(i)
	def test_membership(self,kmers_test: tuple):
		"""
		Method to test if new group of kmers belongs to currently existing group
		or if a new bloom_object needs to be created
		"""
		for i in kmers_test[1]:
			try:
				assert i not in self.bloom
			except AssertionError:
				self.name_list.append(kmers_test[0])
				for j in kmers_test[1]:
					self.bloom.add(j)
				return 1
			return 0

class collapsed_contigs:
	"""
	Object to hold list of bloom_objects
	"""
	def __init__(self):
		self.bloom_list = []
	@staticmethod
	def kmer_assembler(fasta: tuple,kmer: int):
		"""
		Method to chunk out kmers from fasta
		"""
		kmer_list = []
		if len(fasta[1]) < kmer:
			return None
		else:
			for i in range(0,len(fasta[1])-kmer+1):
				kmer_list.append(fasta[1][i:i+kmer].upper())
			return (fasta[0],kmer_list)
	def clump(self,fasta: tuple):
		"""
		Method to iterate through list of bloom filters and add to end if not
		found in list
		"""
		kmers = self.kmer_assembler(fasta,kmer_len)
		if kmers is not None:
			if len(self.bloom_list) < 1:
				self.bloom_list.append(bloom_object(kmers))
			else:
				is_new = 0
				for i in self.bloom_list:
					is_new += i.test_membership(kmers)
				if is_new == 0:
					self.bloom_list.append(bloom_object(kmers))

def main():
	contigs = collapsed_contigs()
	with open(sys.argv[1]) as r:
		j = 0
		for line in r:
			name = line.strip().strip(">")
			dna = next(r).strip()
			contigs.clump((name,dna))
			print([i.name_list for i in contigs.bloom_list if len(i.name_list) > 3])
			print(j)
			j += 1
if __name__ == "__main__":
	main()
