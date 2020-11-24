# function:make a count table out of htseq results
# input: htseq folder; sample name_match; mRNA id
# output: gene and reads count in each sample;sample_info
import subprocess
import sys
import re



#global varible definition
folder="/Users/yiranli/Dropbox/master_assembly/count"
sample_info="/Users/yiranli/Dropbox/master_assembly/sample_info"
output="/Users/yiranli/Dropbox/master_assembly/all_count.tab"
name_match="/Users/yiranli/Dropbox/master_assembly/sample_time_match"
mRNA_id="/Users/yiranli/Dropbox/master_assembly/CryptoDB-43_CparvumIowaII_mRNA.id"
command="ls "+folder+"|grep count"
p=subprocess.Popen(command,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
file_name=p.communicate()[0].decode("utf-8")
file_name=list(filter(None,file_name.split("\n")))
library=["IPEC","Gal","MDBK","organoid_SI","Upenn"]
time_point=["oocyst","_0h","_2h","_24h","_48h","_72h"]
dict_count={}
mRNA_list=[]
name_dict={}
dict_sample_info={}

with open(name_match) as f:
	for line in f:
		l_features=line.strip().split()
		l_time_point=l_features[0]
		l_library=l_features[2]
		l_string=l_features[1]
		name_dict[l_string]=[l_time_point,l_library]

#function definition
def get_sample_info(a_file_name):
	file_name_library=""
	file_name_time=""
	a=0#control each sample only belong to one library
	b=0#control each sample only belong to one time point
	for x_library in library:
		if x_library in a_file_name:
			a=a+1
			file_name_library=x_library
	for y_time_point in time_point:
		if y_time_point in a_file_name:
			b=b+1
			file_name_time=y_time_point
	if a>1 or b>1:
		print("Error:sample {} belongs to multiple time or library:a={},b={}".format(a_file_name,a,b))
		sys.exit()
	if a==0:
		file_name_library=file_library_search(a_file_name)
	if b==0:
		file_name_time=file_time_search(a_file_name)
	return file_name_time, file_name_library


def file_time_search(a_file_name):
	file_time=""
	for i in name_dict.keys():
		if i in a_file_name:
			file_time=name_dict[i][0]
	if file_time=="":
		print("sample {} time point not found".format(a_file_name))
		sys.exit()
	return file_time

def file_library_search(a_file_name):
	file_library=""
	for i in name_dict.keys():
		if i in a_file_name:
			file_library=name_dict[i][1]
	if file_library=="":
		print("sample {} library not found".format(a_file_name))
		sys.exit()
	return file_library



for i_file_name in file_name:
	file_fullpath=folder+"/"+i_file_name
	i_file_name_time,i_file_name_library=get_sample_info(i_file_name)
	if i_file_name_time.startswith("_"):
		i_file_name_time=i_file_name_time[1:]
	i_file_time_library=i_file_name_time+"_"+i_file_name_library
	dict_sample_info[i_file_time_library]=[i_file_name_time,i_file_name_library]
	with open(file_fullpath) as f:
		for line in f:
			if not line.startswith("_"):
				features=line.rstrip().split()
				gene_id=features[0]+"-RA"
				count=features[1]
				if gene_id in dict_count.keys():
					if i_file_time_library in dict_count[gene_id].keys():
						dict_count[gene_id][i_file_time_library].append(count)
					else:
						dict_count[gene_id][i_file_time_library]=[count]
				else:
					dict_count[gene_id]={i_file_time_library:[count]}

with open(mRNA_id) as f:
	for line in f:
		l_name=line.strip()
		mRNA_list.append(l_name)

print_order=["oocyst_IPEC","oocyst_Gal","0h_Gal","0h_Upenn","2h_Gal","2h_MDBK","24h_MDBK","24h_organoid_SI","24h_Upenn","24h_IPEC","48h_MDBK","48h_Upenn","72h_organoid_SI"]
first_line="gene_id"


#print 
with open(output,"w") as write_count:
	with open(sample_info,"w") as write_info:
		write_info.write("time\tgroup\n")
		for i_print_order in print_order:
			print_tmp_a=['{}_s{}'.format(a,b) for a,b in zip([i_print_order]*len(dict_count["cgd8_590-RA"][i_print_order]),[i for i in range(1,len(dict_count["cgd8_590-RA"][i_print_order])+1)])]
			first_line=first_line+"\t"+"\t".join(print_tmp_a)
			for n in range(len(dict_count["cgd8_590-RA"][i_print_order])):
				write_info.write('{}\t{}\n'.format(dict_sample_info[i_print_order][0],dict_sample_info[i_print_order][1]))
		write_count.write("{}\n".format(first_line))
		for gene_id in dict_count.keys():
			print_tmp=""
			for i_print_order in print_order:
				print(gene_id)
				if i_print_order!="72h_organoid_SI":
					print_tmp=print_tmp+"\t".join(dict_count[gene_id][i_print_order])+"\t"
				else:
					print_tmp=print_tmp+"\t".join(dict_count[gene_id][i_print_order])
			if gene_id in mRNA_list:
				write_count.write('mRNA_{}\t{}\n'.format(gene_id,print_tmp))
			else:
				write_count.write('{}\t{}\n'.format(gene_id,print_tmp))








