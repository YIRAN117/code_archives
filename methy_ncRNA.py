#bedtools intersect -a /Users/yiranli/Dropbox/lncRNA_article/submit_ATCC_oct2020.gff -b /Users/yiranli/Dropbox/lncRNA_article/submit/Cparvum_revisited_ncrna_v3.gff -s -wa -wb -f 0.4|awk '{print $9"\t"$18}'|grep -v "Parent"|grep -v "unspecified"|sed 's/=/    /'g|sed 's/;/    /g'|sed 's/(.)/-RA/g'|awk '{print $2"\t"$4}'>lnc_CPATCC_name_match
# gene_dict={}
# name_match={}
# import re
# with open("/Users/yiranli/Dropbox/lncRNA_article/methylation/lnc_CPATCC_name_match") as f:
# 	for line in f:
# 		features=line.split()
# 		name_match[features[0]]=features[1]

# with open("/Users/yiranli/Dropbox/lncRNA_article/submit/Cparvum_revisited_ncrna_v3.gff") as f1:
# 	for line in f1:
# 		if not line.startswith("#"):
# 			features=line.split("\t")
# 			if features[2]=="gene":
# 				name=re.search(r'ID=.*;descri',features[8]).group()
# 				name=name.replace("ID=",'').replace(";descri","")
# 				gene_dict[name]=[features[0],features[3],features[4],features[6]]
# with open("/Users/yiranli/Dropbox/lncRNA_article/submit_ATCC_oct2020.gff") as f2:
# 	for line in f2:
# 		if not line.startswith("#"):
# 			features=line.split("\t")
# 			if features[2]=="gene":
# 				name=re.search(r'ID=.*\(',features[8]).group()
# 				name=name.replace("ID=",'').replace("(","-RA")
# 				gene_dict[name]=[features[0],features[3],features[4],features[6]]
# with open("/Users/yiranli/Dropbox/lncRNA_article/methylation/info") as f3:
# 	for line in f3:
# 		if not line.startswith("Gene_ID"):
# 			features=line.split("\t")
# 			lnc=features[0]
# 			upstream=features[16]
# 			cor=features[17]
# 			try:
# 				if gene_dict[lnc][3]!=gene_dict[upstream][3]:
# 					if gene_dict[lnc][3]=="+":
# 						print('{}\t{}\t{}\t{}A{}\t{}\t{}'.format(gene_dict[lnc][0],int(gene_dict[lnc][1])-200,int(gene_dict[lnc][1]),lnc,name_match[lnc],cor,gene_dict[lnc][3]))
# 					else:
# 						print('{}\t{}\t{}\t{}A{}\t{}\t{}'.format(gene_dict[lnc][0],int(gene_dict[lnc][2]),int(gene_dict[lnc][2])+200,lnc,name_match[lnc],cor,gene_dict[lnc][3]))
# 				else:
# 					if gene_dict[lnc][3]=="+":
# 						print('{}\t{}\t{}\t{}B{}\t{}\t{}'.format(gene_dict[lnc][0],int(gene_dict[lnc][1])-200,int(gene_dict[lnc][1]),lnc,name_match[lnc],cor,gene_dict[lnc][3]))
# 					else:
# 						print('{}\t{}\t{}\t{}B{}\t{}\t{}'.format(gene_dict[lnc][0],int(gene_dict[lnc][2]),int(gene_dict[lnc][2])+200,lnc,name_match[lnc],cor,gene_dict[lnc][3]))					
# 			except:
# 				if gene_dict[lnc][3]=="+":
# 					print('{}\t{}\t{}\t{}C{}\t{}\t{}'.format(gene_dict[lnc][0],int(gene_dict[lnc][1])-200,int(gene_dict[lnc][1]),lnc,name_match[lnc],cor,gene_dict[lnc][3]))
# 				else:
# 					print('{}\t{}\t{}\t{}C{}\t{}\t{}'.format(gene_dict[lnc][0],int(gene_dict[lnc][2]),int(gene_dict[lnc][2])+200,lnc,name_match[lnc],cor,gene_dict[lnc][3]))				


#python /Users/yiranli/Dropbox/scripts/methy_ncRNA.py >lnc_200.bed
#bedtools intersect -a CHH_all.bed -b lnc_200.bed -wa -wb -f 0.9>bed_200.bed
# with open("/Users/yiranli/Dropbox/lncRNA_article/methylation/bed_200.bed") as f4:
# 	methy_dict={}
# 	for line in f4:
# 		features=line.split()
# 		methy=float(features[3])
# 		try:
# 			cor=float(features[8])
# 		except:
# 			cor=0
# 		lnc=features[7]
# 		if lnc not in methy_dict.keys():
# 			methy_dict[lnc]=[methy,cor]
# 		else:
# 			methy_dict[lnc]=[methy+methy_dict[lnc][0],cor]
# 	for i in methy_dict.keys():
#		print(i,methy_dict[i][0],methy_dict[i][1],sep="\t")

#python /Users/yiranli/Dropbox/scripts/methy_ncRNA.py|sort -k1,1nr|sed 's/-RAA/    head_to_tail    /g'|sed 's/-RAB/    head_to_head    /g'|sed 's/-RAC/    NA    /g'

dict1={}
with open("/Users/yiranli/Dropbox/lncRNA_article/methylation/CHH.tab") as f1:
	for line in f1:
		features=line.split()
		dict1[features[0]]=line.rstrip()
dict2={}
with open("/Users/yiranli/Dropbox/lncRNA_article/methylation/CHG.tab") as f2:
	for line in f2:
		features=line.split()
		dict2[features[0]]=features[3]
dict3={}
with open("/Users/yiranli/Dropbox/lncRNA_article/methylation/CHG.tab") as f3:
	for line in f3:
		features=line.split()
		dict3[features[0]]=features[3]

for i in dict1.keys():
	if i in dict2.keys():
		dict1[i]='{}\t{}'.format(dict1[i],dict2[i])
	else:
		dict1[i]='{}\t{}'.format(dict1[i],"NA")
	if i in dict3.keys():
		dict1[i]='{}\t{}'.format(dict1[i],dict3[i])
	else:
		dict1[i]='{}\t{}'.format(dict1[i],"NA")


for i in dict1.keys():
	print(dict1[i])




















	