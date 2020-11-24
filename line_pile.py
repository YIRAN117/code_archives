#function: display gene reative postion (lncRNA overlapped positon with mRNAs)
#input: two gff/bed files that is going to be compared
#output: normalized gene relative postion (from -50%-150% )
import argparse
from pybedtools import BedTool
import re

# parser = argparse.ArgumentParser(description='function: display gene reative postion (lncRNA overlapped positon with mRNAs)')
# parser.add_argument('file1', type=str, help='mRNA gff/bed file')
# parser.add_argument('file2', type=str, help='lncRNA gff/bed file')
# args = parser.parse_args()

lncrna_db=BedTool("/Volumes/Seagate/small_rna/ATCC_genome/CryptoDB-46_CparvumIOWA-ATCC_mRNA.gff")
mrna_db=BedTool("/Volumes/Seagate/small_rna/ATCC_genome/CryptoDB-46_CparvumIOWA-ATCC_mRNA.gff")
antisense=mrna_db.intersect(lncrna_db,wao=True)
a=0
output=[]
print("start\tend\tstrand\tindex")
for entry in antisense:
	if entry[18]!="0":
		a=a+1
		mRNA_id=re.search("CPATCC_\d+-RA",entry[8]).group()
		if entry[6]=="+":
			mRNA_start=entry[3]
			mRNA_end=entry[4]
		else:
			mRNA_start=entry[4]
			mRNA_end=entry[3]
		lncRNA_id=re.search("CPATCC_\d+-RA",entry[17]).group()
		if entry[15]=="+":
			lncRNA_start=entry[12]
			lncRNA_end=entry[13]
		else:
			lncRNA_start=entry[13]
			lncRNA_end=entry[12]
		#(a-b)/x*100
		x=abs(int(mRNA_start) - int(mRNA_end))
		if entry[6]=="+":
			start=(int(lncRNA_start) - int(mRNA_start))/x*100
			end=(int(lncRNA_end) - int(mRNA_start))/x*100
		else:
			start=-1*(int(lncRNA_start) - int(mRNA_start))/x*100
			end=-1*(int(lncRNA_end) - int(mRNA_start))/x*100
			middle=(start+end)/2
		if entry[15]==entry[6]:
			strand=1
		else:
			strand=0
		print('{}\t{}\t{}\t{}'.format(min(start,end),max(start,end),strand,a))
		#print(entry)


		
