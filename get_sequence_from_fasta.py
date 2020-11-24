from Bio import SeqIO
fasta_file=input("where is the fasta file?\n")
para1=input("forward ot reverse extract?\n")
name_file=input("where is the list file?\n")
output=input("where is the output file?\n")
name=open(name_file).readlines()
name=[n.replace("\n","") for n in name]
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
if para1=="forward":
	with open(output,"w") as f:
		for fasta in fasta_sequences:
			nn,ss=fasta.id,fasta.seq
			if nn in name:
				f.write('>{}\n{}\n'.format(nn,ss))
elif para1=="reverse":
	with open(output,"w") as f:
		for fasta in fasta_sequences:
			nn,ss=fasta.id,fasta.seq
			if nn not in name:
				f.write('>{}\n{}\n'.format(nn,ss))


