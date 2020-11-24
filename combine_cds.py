from sys import argv
from Bio import SeqIO
input_fasta=argv[1]
out=argv[2]
t=0
cds_combine=""
fasta_sequences = SeqIO.parse(open(input_fasta),'fasta')
for fasta in fasta_sequences:
	sequence=fasta.seq
	if sequence.startswith("ATG") and sequence.endswith(("TAA", "TAG", "TGA")) and (len(sequence) % 3 ==0):
		t+=1
		cds_combine=cds_combine+sequence

with open(out,"w") as f2:
	f2.write(">combine.cds\n{}".format(cds_combine))

print("complete cds: {}".format(t))


