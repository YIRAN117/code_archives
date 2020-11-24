import re
pro_db={}
a=0
b=0
c=0
d=0
e=0
f=0
with open("/home/yl22225/tcru/family/gene_family/Brazil/protein/brazil_protein.list") as f1:
	for line in f1:
		pro_id=line.rstrip()
		pro_db[pro_id]="hypothetical protein"
		a=a+1
	print("total number of protein is {}".format(a))

with open("/home/yl22225/tcru/family/gene_family/Brazil/blast_t/out_Non-Esmeraldo-like.name") as f2:
	for line in f2:
		features=line.rstrip().split("\t")
		pro_id=features[0]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name and pro_id in pro_db.keys():
			b=b+1
			pro_db[pro_id]=name
	print("total number of naming by Non-Esmeraldo is {}".format(b))
with open("/home/yl22225/tcru/family/gene_family/Brazil/blast_t/out_TcruziCLBrener.name") as f3:
	for line in f3:
		features=line.rstrip().split("\t")
		pro_id=features[0]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name and pro_id in pro_db.keys() and "hypothetical" in pro_db[pro_id]:
			c=c+1
			pro_db[pro_id]=name
	print("total number of naming by CL_Brener is {}".format(c))

with open("/home/yl22225/tcru/family/gene_family/Brazil/blast_t/out_TbruceiTREU927.name") as f4:
	for line in f4:
		features=line.rstrip().split("\t")
		pro_id=features[0]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name and pro_id in pro_db.keys() and "hypothetical" in pro_db[pro_id]:
			d=d+1
			pro_db[pro_id]=name
	print("total number of naming by brucei is {}".format(d))

with open("/home/yl22225/tcru/family/gene_family/Brazil/blast_t/out_LmajorFriedlin.name") as f5:
	for line in f5:
		features=line.rstrip().split("\t")
		pro_id=features[0]
		name=features[2].replace("gene_product=","")
		if "hypothetical" not in name and pro_id in pro_db.keys() and "hypothetical" in pro_db[pro_id]:
			e=e+1
			pro_db[pro_id]=name
	print("total number of naming by Leishmania is {}".format(e))


for key,value in pro_db.items():
	print("{}\tproduct\t{}".format(key,value))




