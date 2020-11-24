import re
name_list=[]
with open("/Users/yiranli/Downloads/transport_id") as f1:
	for line in f1:
		name_list.append(line.rstrip())

name_dict={}
with open("/Users/yiranli/Dropbox/master_assembly/ATCC/Cryptodb/CryptoDB-46_CparvumIOWA-ATCC_mRNA.gff") as f2:
	for line in f2:
		name=re.search("CPATCC_\d+",line).group()
		description=re.search("description=.*$",line).group().replace("description=","")
		name_dict[name]=description

for i in name_list:
	print(name_dict[i])

