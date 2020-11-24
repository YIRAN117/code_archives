#function: make long fasta file into small parts allowed in online tools
#input: one fasta file (one file or multiple)
#output: multiple fasta in the same folder
import argparse
import re
parser = argparse.ArgumentParser(description='make long fasta file into small parts allowed in online tools')
parser.add_argument('input', type=str, 
                    help='one fasta file (one file or multiple)')
parser.add_argument('size', type=int, 
                    help='split every n string')
args = parser.parse_args()
n=args.size
with open(args.input) as f:
	all_content=f.read()
	all_content=re.sub('>.*',"",all_content)
	all_content=re.sub('\n',"",all_content)
	out=[all_content[i:i+n] for i in range(0,len(all_content),n)]
	for i,j in enumerate(out):
		name=f'{i}.fa'
		with open(name,"w") as write_to:
			write_to.write('>{}\n{}\n'.format(i,j))