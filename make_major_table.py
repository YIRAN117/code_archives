# function: integrate information and make the major table "ID oocyst DE? h0 DE? h24 DE? antisense antisense_cor antisense_p.. "
# input: DE folder, vst_table,
# output: the major table



g_input_DE_folder="/Users/yiranli/Dropbox/lncRNA_article/data/DE"
g_input_vst_table="/Users/yiranli/Dropbox/lncRNA_article/data/all_count_vsd_table.tab"
g_output_major_table="/Users/yiranli/Dropbox/lncRNA_article/data/major_table.tab"
import os
import subprocess
from statistics import mean,stdev

#make DE diction
g_DE_dict={}
os.chdir(g_input_DE_folder)
g_p=subprocess.Popen(["ls", "."],stdout=subprocess.PIPE)
g_DE_timename,error=g_p.communicate()

if error:
	raise Exception("DE file name extraction error\n")


g_DE_timename_list=g_DE_timename.decode().split()
for i_timename in g_DE_timename_list:
	l_DE_file=g_input_DE_folder+"/"+i_timename
	l_DE_file.replace("//","/")#in case of redundant folder sign
	with open(l_DE_file) as handle:
		#header=handle.readline()
		for line in handle:
			line=line.rstrip().split()
			name=line[0]
			value=line[1]
			if name in g_DE_dict.keys():
				g_DE_dict[name][i_timename]=value
			else:
				g_DE_dict[name]=dict(zip(g_DE_timename,[0]*len(g_DE_timename_list)))



#make expression diction
g_expression_dict={}
with open(g_input_vst_table) as handle:
	header=handle.readline()
	for line in handle:
		features=line.rstrip().split()
		gene=features[0]
		values=[float(x) for x in features[1:]]
		oocyst=mean(values[0:2])
		oocyst_std=stdev(values[0:2])
		h0=mean(values[2:9])
		h0_std=stdev(values[2:9])
		h2=mean(values[7:19])
		h2_std=stdev(values[7:19])
		h24=mean(values[19:26])
		h24_std=stdev(values[19:26])
		h48=mean(values[26:33])
		h48_std=stdev(values[26:33])
		if gene not in g_expression_dict.keys():
			g_expression_dict[gene]=[oocyst,oocyst_std,h0,h0_std,h2,h2_std,h24,h24_std,h48,h48_std]
		else:
			raise Exception("duplicate gene id in vst table")

#make cor_relation diction
from up_down_anti_cor import dict_upstream,dict_downstream,dict_antisense
from cor_distance import g_cor_dict,g_pvalue_dict

g_cor_relation_dict={}

for gene in g_expression_dict.keys():
	g_cor_relation_dict[gene]={"antisense_gene":"NA","antisense_cor":"NA","antisense_p":"NA","upstream_gene":"NA","upstream_cor":"NA","upstream_p":"NA","downstream_gene":"NA","downstream_cor":"NA","downstream_p":"NA"}


for gene in g_expression_dict.keys():
	if not gene.startswith("mRNA"):					
		try:			
			upstream_gene="mRNA_"+dict_upstream[gene]
			g_cor_relation_dict[gene]["upstream_gene"]=upstream_gene
			g_cor_relation_dict[gene]["upstream_cor"]=g_cor_dict[(gene,upstream_gene)]
			g_cor_relation_dict[gene]["upstream_p"]=g_pvalue_dict[(gene,upstream_gene)]
			g_cor_relation_dict[upstream_gene]["upstream_gene"]=gene
			g_cor_relation_dict[upstream_gene]["upstream_cor"]=g_cor_dict[(gene,upstream_gene)]
			g_cor_relation_dict[upstream_gene]["upstream_p"]=g_pvalue_dict[(gene,upstream_gene)]
		except:
			pass
		try:			
			downstream_gene="mRNA_"+dict_downstream[gene]
			g_cor_relation_dict[gene]["downstream_gene"]=downstream_gene
			g_cor_relation_dict[gene]["downstream_cor"]=g_cor_dict[(gene,downstream_gene)]
			g_cor_relation_dict[gene]["downstream_p"]=g_pvalue_dict[(gene,downstream_gene)]
			g_cor_relation_dict[downstream_gene]["downstream_gene"]=gene
			g_cor_relation_dict[downstream_gene]["downstream_cor"]=g_cor_dict[(gene,downstream_gene)]
			g_cor_relation_dict[downstream_gene]["downstream_p"]=g_pvalue_dict[(gene,downstream_gene)]
		except:
			pass
		try:			
			antisense_gene="mRNA_"+dict_antisense[gene]
			g_cor_relation_dict[gene]["antisense_gene"]=antisense_gene
			g_cor_relation_dict[gene]["antisense_cor"]=g_cor_dict[(gene,antisense_gene)]
			g_cor_relation_dict[gene]["antisense_p"]=g_pvalue_dict[(gene,antisense_gene)]
			g_cor_relation_dict[antisense_gene]["antisense_gene"]=gene
			g_cor_relation_dict[antisense_gene]["antisense_cor"]=g_cor_dict[(gene,antisense_gene)]
			g_cor_relation_dict[antisense_gene]["antisense_p"]=g_pvalue_dict[(gene,antisense_gene)]
		except:
			pass


#write out the table
with open(g_output_major_table,"w") as write_to:
	write_to.write("gene\toocyst\toocyst_std\th0\th0_std\th2\th2_std\th24\th24_std\th48\th48_std\th0_h24_DE\th24_h48_DE\tantisense_gene\tantisense_cor\tantisense_p\tupstream_gene\tupstream_cor\tupstream_p\tdownstream_gene\tdownstream_cor\tdownstream_p\n")
	for gene in g_expression_dict.keys():
		gene=gene.replace("Cgd","cgd")
		try:
			write_to.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene,"\t".join(map(str,g_expression_dict[gene])),g_DE_dict[gene]["DE_h0_h24"],g_DE_dict[gene]["DE_h24_h48"],g_cor_relation_dict[gene]["antisense_gene"],g_cor_relation_dict[gene]["antisense_cor"],g_cor_relation_dict[gene]["antisense_p"],g_cor_relation_dict[gene]["upstream_gene"],g_cor_relation_dict[gene]["upstream_cor"],g_cor_relation_dict[gene]["upstream_p"],g_cor_relation_dict[gene]["downstream_gene"],g_cor_relation_dict[gene]["downstream_cor"],g_cor_relation_dict[gene]["downstream_p"]))
		except:
			pass













