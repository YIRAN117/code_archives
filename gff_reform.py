import re
with open("/Users/yiranli/Desktop/intron/Annotations.gff3") as f:
	with open("/Users/yiranli/Desktop/intron/reformed.gff3","w") as write_to:
	#head=[next(f) for x in range(10)]
		dict={}
		a=0
		for rline in f:
			a+=1
			print(a)
			if not rline.startswith("#"):
				line=rline.rstrip()+";"
				print(line)
				feature=line.rstrip().split("\t")
				name=re.search(r"Name=.+?;",feature[8]).group()
				Id=re.search(r"ID=.+?;",feature[8]).group()
				try:
					parent=re.search(r"Parent=.+?;",feature[8]).group()
				except:
					print("no parent")
				key=Id
				val=name
				dict[key]=val
				if feature[2]=="gene":
					last_col=name.replace("Name","ID").replace(";","")
				elif feature[2]=="mRNA":
					parent=dict[parent.replace("Parent","ID")]
					last_col=name.replace("Name","ID")+parent.replace("Name","Parent").replace(";","")
				elif feature[2]=="tRNA":
					parent=dict[parent.replace("Parent","ID")]
					last_col=name.replace("Name","ID").replace(";","-RA;")+parent.replace("Name","Parent").replace(";","")
				elif feature[2]=="exon" or "CDS":
					parent=dict[parent.replace("Parent=","ID=")]
					last_col=parent.replace("Name=","Parent=")
				write_to.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(feature[0],feature[1],feature[2],feature[3],feature[4],feature[5],feature[6],feature[7],last_col))


			
