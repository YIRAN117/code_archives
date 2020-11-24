#starts from only gene features
#usage python3.6 exonerategff_reform.py input.gff
import re
import sys
file=sys.argv[1]
out=file+".gff3"
with open(file) as f:
	with open(out,"w") as write_to:
		a=0
		for rline in f:
			a+=1
			line=rline.rstrip()
			feature=line.rstrip().split("\t")
			name=re.search(r"sequence.*;",feature[8]).group().replace("sequence","").replace(" ","").replace(";","")+"-"+str(a)
			ID_gene=name+":gene"
			ID_mRNA=name+":mRNA"
			ID_exon=name+":exon"
			ID_CDS=name+":CDS"
			gene='{}\t{}\tgene\t{}\t{}\t{}\t{}\t{}\tID={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_gene)
			mRNA='{}\t{}\tmRNA\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_mRNA,ID_gene)
			exon='{}\t{}\texon\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_exon,ID_mRNA)
			CDS='{}\t{}\tCDS\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_CDS,ID_mRNA)
			write_to.write('{}{}{}{}'.format(gene,mRNA,exon,CDS))