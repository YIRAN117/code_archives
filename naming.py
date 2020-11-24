import re
pro_db={}
a=0
b=0
c=0
d=0
e=0
f=0
with open("/Users/yiranli/Desktop/brazil/last_2019/blast_t/protein_list") as f1:
	for line in f1:
		pro_id=line.rstrip()[0:36]
		pro_db[pro_id]="hypothetical protein"
		a=a+1
	print("total number of protein is {}".format(a))

with open("/Users/yiranli/Desktop/brazil/last_2019/blast_t/1e-10/Trypanosoma_cruzi_CL_Brener_Non-Esmeraldo-like.name") as f2:
	for line in f2:
		features=line.rstrip().split("\t")
		pro_id=features[0].replace("protein_","")[0:36]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name:
			b=b+1
			pro_db[pro_id]=name
	print("total number of naming by Non-Esmeraldo is {}".format(b))
with open("/Users/yiranli/Desktop/brazil/last_2019/blast_t/1e-10/Trypanosoma_cruzi_strain_CL_Brener.name") as f3:
	for line in f3:
		features=line.rstrip().split("\t")
		pro_id=features[0].replace("protein_","")[0:36]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name and "hypothetical" in pro_db[pro_id]:
			c=c+1
			pro_db[pro_id]=name
	print("total number of naming by CL_Brener is {}".format(c))

with open("/Users/yiranli/Desktop/brazil/last_2019/blast_t/1e-10/Trypanosoma_brucei.name") as f4:
	for line in f4:
		features=line.rstrip().split("\t")
		pro_id=features[0].replace("protein_","")[0:36]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name and "hypothetical" in pro_db[pro_id]:
			d=d+1
			pro_db[pro_id]=name
	print("total number of naming by brucei is {}".format(d))

with open("/Users/yiranli/Desktop/brazil/last_2019/blast_t/1e-10/Leishmania_major.name") as f5:
	for line in f5:
		features=line.rstrip().split("\t")
		pro_id=features[0].replace("protein_","")[0:36]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name and "hypothetical" in pro_db[pro_id]:
			e=e+1
			pro_db[pro_id]=name
	print("total number of naming by Leishmania is {}".format(e))


with open("/Users/yiranli/Desktop/brazil/last_2019/blast_t/1e-10/nr.name") as f6:
	for line in f6:
		features=line.rstrip().split("\t")
		pro_id=features[0].replace("protein_","")[0:36]
		name=re.sub(r"\[.*\]","",features[1])
		if "hypothetical" not in name and "hypothetical" in pro_db[pro_id]:
			f=f+1
			pro_db[pro_id]=name
	print("total number of naming by nr is {}".format(f))

for key,value in pro_db.items():
	print("{}\tproduct\t{}".format(key,value))




