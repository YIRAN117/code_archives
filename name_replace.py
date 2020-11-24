import re
with open("/Users/yiranli/Dropbox/lncRNA_article/old_new_names") as f1:
	name_dict={lines.rstrip().split("\t")[0]:lines.rstrip().split("\t")[1] for lines in f1}
#print(name_dict)
with open("/Users/yiranli/Dropbox/lncRNA_article/data/mRNA_cluster.tab") as f2:
	for lines in f2:
		name=re.findall(r'cgd\d_\d+-RA',lines)
		for i in name:
			#print(i)
			try:
				lines=lines.replace(i,name_dict[i]).replace("\n","")
			except:
				pass
		print(lines)




# import re
# with open("/Users/yiranli/Dropbox/lncRNA_article/lncRNA_name_match") as f1:
# 	name_dict={lines.rstrip().split("\t")[0]:lines.rstrip().split("\t")[1] for lines in f1}
# print(name_dict)
# with open("/Users/yiranli/Dropbox/lncRNA_article/conservation/common167_lncRNA.bed") as f2:
# 	for lines in f2:
# 		name=lines.split("\t")
# 		name=name[8].replace("Parent=","").replace("\n","")
# 		try:
# 			lines=lines.replace(name,name_dict[name]).replace("\n","")
# 		except:
# 			pass
# 		print(lines)