import re
name_dict={}
with open("/home/yl22225/tcru/family/gene_family/Brazil/naming/pseudogene/name_result") as f1:
	for lines in f1:
		name_dict[lines.split("\t")[0]]=lines.split("\t")[2].rstrip()
# for key,values in name_dict.items():
# 	print(key,values)

with open("/home/yl22225/tcru/family/gene_family/Brazil/naming/pseudogene/brazil_pseudogene.gff3") as f2:
	for lines in f2:
		features=lines.split("\t")
		if features[2]=="pseudogene":
			name=features[8].split(":")[0].replace("ID=","")
			lines=lines.rstrip()+";Name="+re.sub(r"^ ","",name_dict[name])
			print(lines)
		else:
			print(lines.rstrip())