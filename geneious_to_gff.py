import re
import sys
file=sys.argv[1]
out=file+".gff3"
with open(file) as f:
	with open(out,"w") as write_to:
		a=0
		for rline in f:
			line=rline.rstrip()
			feature=line.rstrip().split("\t")
			if feature[2]=="CDS":
				print(feature)
				try:
					name=re.search(r"ID=.*;",feature[8]).group().replace("ID=","").replace(" ","").replace(";","")
				except:
					print("error")
					name="wei_curated"+str(a)
					a=a+1
				ID_gene=name+":gene"
				ID_mRNA=name+":mRNA"
				ID_exon=name+":exon"
				ID_CDS=name+":CDS"
				gene='{}\t{}\tgene\t{}\t{}\t{}\t{}\t{}\tID={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_gene)
				mRNA='{}\t{}\tmRNA\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_mRNA,ID_gene)
				exon='{}\t{}\texon\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_exon,ID_mRNA)
				CDS='{}\t{}\tCDS\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_CDS,ID_mRNA)
				write_to.write('{}{}{}{}'.format(gene,mRNA,exon,CDS))
			elif feature[2]=="pseudogene":
				try:
					name=re.search(r"ID=.*",feature[8]).group().replace("ID=","").replace(" ","").replace(";","")
				except:
					name="wei_curated"+str(a)
					a=a+1
				ID_gene=name+":pseudogene"
				ID_mRNA=name+":transcript"
				ID_exon=name+":exon"
				pseudogene='{}\t{}\tpseudogene\t{}\t{}\t{}\t{}\t{}\tID={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_gene)
				transcript='{}\t{}\ttranscript\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_mRNA,ID_gene)
				exon='{}\t{}\texon\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};\n'.format(feature[0],feature[1],feature[3],feature[4],feature[5],feature[6],feature[7],ID_exon,ID_mRNA)
				write_to.write('{}{}{}'.format(pseudogene,transcript,exon))