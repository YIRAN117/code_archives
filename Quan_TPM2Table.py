# function:make a TPM table out of stringtie -e results
# input: stringtie -e folder; sample name_match; mRNA id
# outpur: gene and TPM level in each sample 
import subprocess
import sys
import re


#global varible definition
folder="/Users/yiranli/Dropbox/master_assembly/quan"
name_match="/Users/yiranli/Dropbox/master_assembly/sample_time_match"
mRNA_id="/Users/yiranli/Dropbox/master_assembly/CryptoDB-43_CparvumIowaII_mRNA.id"
command="ls "+folder+"|grep abundance"
p=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
file_name=p.communicate()[0].decode("utf-8")
file_name=list(filter(None,file_name.split("\n")))
library=["IPEC","Gal","MDBK","organoid_SI","upenn"]
time_point=["oocyst","_0h","_2h","_24h","_48h","_72h"]
dict_TPM={}
mRNA_list=[]

#function definition
def get_sample_info(file_name):
	file_name_library=""
	file_name_time=""
	a=0#control each sample only belong to one library
	b=0#control each sample only belong to one time point
	for x_library in library:
		if x_library in file_name:
			a=a+1
			file_name_library=x_library
			for y_time_point in time_point:
				if y_time_point in file_name:
					b=b+1
					file_name_time=y_time_point
	if a>1 or b>1:
		print("Error:sample {} belongs to multiple time or library".format(i_file_name))
		sys.exit()
	if a==0:
		print("Error:sample {} library not found".format(i_file_name))
		sys.exit()
	if b==0:
		file_name_time=file_name_search(file_name,name_match)

	return file_name_time, file_name_library


def file_name_search(file_name, name_match):
	name_dict={}
	file_time=""
	with open(name_match) as f:
		for line in f:
			features=line.strip().split()
			time_point=features[0]
			string=features[1]
			name_dict[string]=time_point
	for i in name_dict.keys():
		if i in file_name:
			file_time=name_dict[i]
	if file_time=="":
		print("sample {} time point not found",format(file_name))
		sys.exit()
	return file_time


for i_file_name in file_name:
	file_fullpath=folder+"/"+i_file_name
	i_file_name_time,i_file_name_library=get_sample_info(i_file_name)
	i_file_time_library=i_file_name_time.replace("_","",1)+"_"+i_file_name_library
	with open(file_fullpath) as f:
		next(f)
		for line in f:
			if not line.startswith("STRG"):
				features=line.rstrip().split("\t")
				gene_id=features[0]+"-RA"
				TPM=features[8]
				if gene_id in dict_TPM.keys():
					if i_file_time_library in dict_TPM[gene_id].keys():
						dict_TPM[gene_id][i_file_time_library].append(TPM)
					else:
						dict_TPM[gene_id][i_file_time_library]=[TPM]
				else:
					dict_TPM[gene_id]={i_file_time_library:[TPM]}

with open(mRNA_id) as f:
	for line in f:
		name=line.strip()
		mRNA_list.append(name)

print_order=["oocyst_IPEC","oocyst_Gal","0h_Gal","0h_upenn","2h_Gal","2h_MDBK","24h_MDBK","24h_organoid_SI","24h_upenn","24h_IPEC","48h_MDBK","48h_upenn","72h_organoid_SI"]
first_line="gene_id"
for i in print_order:
	first_line=first_line+"\t"+"\t".join([i]*len(dict_TPM["cgd8_590-RA"][i]))
print(first_line)
for gene_id in dict_TPM.keys():
	print_tmp=""
	for i_print_order in print_order:
		if i_print_order!="72h_organoid_SI":
			print_tmp=print_tmp+"\t".join(dict_TPM[gene_id][i_print_order])+"\t"
		else:
			print_tmp=print_tmp+"\t".join(dict_TPM[gene_id][i_print_order])
	if gene_id in mRNA_list:
		print('mRNA_{}\t{}'.format(gene_id,print_tmp))
	else:
		print('{}\t{}'.format(gene_id,print_tmp))






