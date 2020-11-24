import re
anno_product={}
anno_Dbxref={}
anno_Ontology_term={}
with open("/home/yl22225/tcru/family/gene_family/Brazil/naming/protein/genome.annotations") as f1:
	for line in f1:
		features=line.rstrip().split("\t")
		if len(features)!=3:
			pass
		else:
			pro_id=features[0]
			source=features[1]
			info=features[2]				
			if source=="product":
				anno_product[pro_id]=info
			elif source=="Dbxref" and "GO:" in info:
				if pro_id not in anno_Ontology_term.keys():
					anno_Ontology_term[pro_id]=info.replace("|GO:",",GO:")
				else:
					anno_Ontology_term[pro_id]=anno_Ontology_term[pro_id]+","+info.replace("|GO:",",GO:")
			elif source=="Dbxref" and "GO:" not in info:
				if pro_id not in anno_Dbxref.keys():
					anno_Dbxref[pro_id]=info
				else:
					anno_Dbxref[pro_id]=anno_Dbxref[pro_id]+","+info

with open("/home/yl22225/tcru/family/gene_family/Brazil/naming/protein/Brazil_all_protein.gff3") as f2:
	for line in f2:
		features=line.rstrip().split("\t")
		if len(features)!=9:
			print("headers:{}".format(line))
		else:
			if features[2]=="gene":
				pro_id=features[8].split(":")[0].replace("ID=","")
				try:
					features[8]=features[8]+";Name="+anno_product[pro_id]
				except:
					pass
			elif features[2]=="mRNA":
				pro_id=features[8].split(":")[0].replace("ID=","")
				try:
					features[8]=features[8]+";Ontology_term="+anno_Ontology_term[pro_id]
				except:
					pass
				try:
					features[8]=features[8]+";Dbxref="+anno_Dbxref[pro_id]
				except:
					pass				
			else:
				features[8]=features[8]
			print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(features[0],features[1],features[2],features[3],features[4],features[5],features[6],features[7],features[8]))
			



	