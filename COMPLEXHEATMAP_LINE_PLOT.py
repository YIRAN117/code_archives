import argparse

parser = argparse.ArgumentParser(description='function: calculate atgc percent')
parser.add_argument('file', type=str, help='input seekr file')
args = parser.parse_args()

kmer=[]
kmer_dict={}
x=0
index=0

print("A\tT\tG\tC")
with open(args.file) as f1:
	header=f1.readline()
	kmer=list(filter(None,header.split(",")))
	x=len(kmer[0])
	for i in kmer:
		index+=1
		A=i.count("A")/x*100
		T=i.count("T")/x*100
		C=i.count("C")/x*100
		G=i.count("G")/x*100		
		print('{}\t{}\t{}\t{}'.format(A,T,G,C))