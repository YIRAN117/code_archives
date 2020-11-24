import os
g_basepath="/scratch/yl22225/small_rna/rnaz/rnaz/result/"
entries = list(os.scandir(g_basepath))
file_name=[entry.name for entry in entries if entry.name.endswith(".dat")]
for i in file_name:
	with open(os.path.join(g_basepath,i)) as f:
		strand=""
		p=0
		for line in f:
			if line.startswith('window'):
				#print(i,line)
				features=line.rstrip().split("\t")								
				if float(features[17]) >p:
					p=float(features[17])
					strand=features[5]
		if strand in ["-","+"]:
			print(i.replace(".dat",""),p,strand,sep="\t")
		else:
			print(i.replace(".dat",""),p,".",sep="\t")

