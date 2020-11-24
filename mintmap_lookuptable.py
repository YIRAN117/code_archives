from Bio.Seq import *
import re
genome_seq=""
with open("/Users/yiranli/Downloads/genome.fasta") as f1:
	for line in f1:
		if not line.startswith(">"):
			genome_seq=genome_seq+line.strip()
genome_seq_rc=reverse_complement(genome_seq)

with open("/Users/yiranli/Downloads/cpar_tRNA_cryptodb.fa") as f2:
	for line in f2:
		# if not line.startswith(">"):
		# 	print(line.strip())
		# else:
		# 	i_seq=line.strip()
		# 	print(reverse_complement(i_seq))
		if not line.startswith(">"):
			i_seq=line.strip()
			i_seq=reverse_complement(i_seq)
			for j in range(0,len(i_seq)-16):
				for i in range(j+16,len(i_seq)+1):
					count=len(re.findall(i_seq[j:i],genome_seq))+len(re.findall(i_seq[0:i],genome_seq_rc))
					if count==1:
						print('{}\tY'.format(i_seq[j:i]))
					else:
						print('{}\tN'.format(i_seq[j:i]))