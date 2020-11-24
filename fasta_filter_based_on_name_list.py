#usage:python3.6 script.py mRNA.fa name_list
#output:mRNA.fa_trans.out

from Bio import SeqIO
import sys
fasta_file=sys.argv[1]
name_file=sys.argv[2]
name_list=[]
out_file=fasta_file+"_complete.fa"

with open(name_file) as f0:
	for i in f0:
		name_list.append(i.rstrip())
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

with open(out_file,"w") as f1:
	for fasta in fasta_sequences:
		header=fasta.id
		sequence=fasta.seq
		if not any(s in header for s in name_list):
			f1.write(">{}\n{}\n".format(header,sequence))

