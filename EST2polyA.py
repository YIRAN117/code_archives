import re
from Bio import SeqIO
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

g_input_EST="/Users/yiranli/Desktop/EST/EST.fa"
g_EST_polya="/Users/yiranli/Desktop/EST/EST_polya.fa"

with open(g_EST_polya,"w") as f0:
	with open(g_input_EST) as f1:
		for record in SeqIO.parse(f1,"fasta"):
			seq=str(record.seq)
			seq_rc=Seq(seq).reverse_complement()
			seq_=""
			if re.search('A{10,}[TGC]{0,3}$',str(seq)):
				seq_=re.sub('A{10,}$','',str(seq))
			elif re.search('A{10,}$',str(seq_rc)):
				seq_=re.sub('A{10,}$','',str(seq_rc))
			if len(seq_)>20:
				f0.write('>{}\n{}\n'.format(record.id,seq_))

os.system('blastn -num_threads 4 -db /Users/yiranli/Dropbox/master_assembly/lncRNA_identification/lncRNA_cryptoDB_merged.fa -query /Users/yiranli/Desktop/EST/EST_polya.fa -strand "plus"   -evalue 1e-5 -outfmt 5 >/Users/yiranli/Desktop/EST/EST_blast.xml')
os.system('python /Users/yiranli/Dropbox/scripts/blastxml_to_tabular.py -c ext /Users/yiranli/Desktop/EST/EST_blast.xml| sort -k1,1 -k12,12nr -k11,11n|sort -u -k1,1 >/Users/yiranli/Desktop/EST/EST_blast.tab')

with open("/Users/yiranli/Desktop/EST/EST_blast.tab") as f2:
	for line in f2:
		features=line.rstrip().split()
		send=int(features[9])#End of alignment in subject
		slen=int(features[23])#Subject sequence length
		sseqid=features[1]
		qseqid=features[0]
		if send > (slen -15):
			print('{}\t{}\n'.format(sseqid,qseqid))



