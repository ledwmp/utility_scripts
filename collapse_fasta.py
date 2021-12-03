from bloom_filter import BloomFilter
import sys
import os

kmer_len = 25

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
				assert i not in self.bloom #check if not in bloom
			except AssertionError: #if might be in bloom
				#for j in kmers_test[1]:
				#	self.bloom.add(j) #add all kmers to bloom
				return 1
			return 0
	def offload_kmers(self,kmers_test: tuple):
		self.name_list.append(kmers_test[0])
		for i in kmers_test[1]:
			self.bloom.add(i)

class collapsed_contigs:
	"""
	Object to hold and update list of bloom_objects
	"""
	def __init__(self,kmer_len: int,max_elements: int):
		self.kmer_len = kmer_len
		self.bloom_list = []
		self.bloom_filter = BloomFilter(max_elements=max_elements,error_rate=1./self.kmer_len)
		self.bloom_elements = 0
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
		found in list. Also filters through main bloom filter to avoid iterating
		through full list if all kmers haven't been found.
		"""
		kmers = self.kmer_assembler(fasta,self.kmer_len)
		if kmers is not None:
			if len(self.bloom_list) < 1:
				self.bloom_list.append(bloom_object(kmers))
				for i in kmers[1]:
					self.bloom_elements += 1
					self.bloom_filter.add(i)
			else:
				seen = 0
				offload_indices = set()
				for i in kmers[1]:
					try:
						assert i not in self.bloom_filter #when not in large bloom
					except AssertionError: #might be in large bloom
						seen = 1
						for j,k in enumerate(self.bloom_list):
							#this offloads all kmers into bloom_filter, forcing the
							#next iteration to catch and reiterate again. Instead
							#store indices of catches and then offload into all
							#indices after all kmers in fasta have been iterated
							#through
							if j not in offload_indices:
								found = k.test_membership(kmers)
								if found == 1:
									offload_indices.add(j)
						#if is_new == 0:
							#self.bloom_list.append(bloom_object(kmers))
						#same here
						#for j in kmers[1]:
						#	self.bloom_filter.add(j)
						#break
				if seen == 0: #hasn't been seen
					self.bloom_list.append(bloom_object(kmers))
				else:
					for i in offload_indices:
						self.bloom_list[i].offload_kmers(kmers)
				for i in kmers[1]:
					self.bloom_filter.add(i)
					self.bloom_elements += 1




def main():
	byte_size = os.path.getsize(sys.argv[1])
	contigs = collapsed_contigs(kmer_len,byte_size*2)
	with open(sys.argv[1]) as r:
		j = 0
		for line in r:
			name = line.strip().strip(">")
			dna = next(r).strip()
			contigs.clump((name,dna))
			print([(i.name_list[0].split(",")[0],len(i.name_list)) for i in contigs.bloom_list if len(i.name_list) > 1])
			print(j,contigs.bloom_elements,byte_size/2)
			j += 1
if __name__ == "__main__":
	main()
