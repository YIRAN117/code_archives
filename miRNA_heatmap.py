# function:show mirna prediction results with heatmaps
# input: mirna prediction tool results
# output: tables and heatmaps

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from collections import Counter

#input and output files
mirna_file="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/targetscan_70/miR_Family_Info_human.fa"
lncrna_file="/Users/yiranli/Dropbox/master_assembly/list_clean_manually_checked/lncRNA_checked.fa"
targetscan_result="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/targetscan_70/clean_lncRNA_human_mirna_targetscan_result"
targetscan_result_out="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/targetscan_70/clean_lncRNA_human_mirna_targetscan_result.tab"
targetscan_result_out_top100="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/targetscan_70/clean_lncRNA_human_mirna_targetscan_result_top100.tab"
miranda_result="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/miranda/miranda_out"
miranda_result_tab="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/miranda/miranda_out.tab"
miranda_result_tab_top100="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/miranda/miranda_out_top100.tab"
pita_result="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/pita/cleanLncrna_human_pita_results_targets.tab"
pita_result_tab="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/pita/cleanLncrna_human_pita_results_targets.tab.tab"
pita_result_tab_top100="/Users/yiranli/Dropbox/master_assembly/miRNA_binding/pita/cleanLncrna_human_pita_results_targets_top100.tab.tab"

#1)read in all interested mirna and lncrna id
def make_name_list(fasta_file):
	out_list=[]
	with open(fasta_file) as handle:
		for line in handle:
			if line.startswith(">"):
				name=line.rstrip().replace(">","")
				name=re.sub('\(.\)',"",name)
				out_list.append(name)
	return(out_list)

mirna_list=make_name_list(mirna_file)
lnc_list=make_name_list(lncrna_file)

#2)make dict: {(lnc,mirna):score}
g_targetscan_out_info={}
with open(targetscan_result) as f1:
	header=f1.readline()
	for line in f1:
		features=line.split()
		lnc=features[0]
		mirna=features[1]
		nmer=features[8]
		key=(lnc,mirna)
		if "6mer" not in nmer:
			if key not in g_targetscan_out_info.keys():
				g_targetscan_out_info[key]=1
			else:
				g_targetscan_out_info[key]+=1

g_miranda_out_info={}
with open(miranda_result) as f2:
	for line in f2:
		if line.startswith(">>"):
			line=line.replace(">>","")
			features=line.split()
			mirna=features[0]
			lnc=features[1]
			score=float(features[2])*float(features[3])*-1
			key=(lnc,mirna)
			g_miranda_out_info[key]=score

g_pita_out_info={}
with open(pita_result) as f3:
	header=f3.readline()
	for line in f3:
		features=line.split()
		mirna=features[1]
		lnc=features[0]
		score=float(features[3])*-1
		key=(lnc,mirna)
		g_pita_out_info[key]=score


#3)make output tables
with open(targetscan_result_out,"w") as write_to:
	write_to.write("{}\n".format("\t".join(lnc_list)))
	for i in mirna_list:
		alist=[]
		for j in lnc_list:
			try:
				alist.append(g_targetscan_out_info[(j,i)])
			except:
				alist.append(0)
		write_to.write("{}\t{}\n".format(i,"\t".join(map(str,alist))))

with open(miranda_result_tab,"w") as write_to2:
	write_to2.write("{}\n".format("\t".join(lnc_list)))
	for i in mirna_list:
		alist=[]
		for j in lnc_list:
			try:
				alist.append(g_miranda_out_info[(j,i)])
			except:
				alist.append(0)
		write_to2.write("{}\t{}\n".format(i,"\t".join(map(str,alist))))

with open(pita_result_tab,"w") as write_to3:
	write_to3.write("{}\n".format("\t".join(lnc_list)))
	for i in mirna_list:
		alist=[]
		for j in lnc_list:
			try:
				alist.append(g_miranda_out_info[(j,i)])
			except:
				alist.append(0)
		write_to3.write("{}\t{}\n".format(i,"\t".join(map(str,alist))))

#4)show heatmaps
#targetscan
targetscan_df=pd.read_csv(targetscan_result_out,index_col=0,sep='\t')
ax = sns.heatmap(
    targetscan_df, 
    vmin=0, vmax=1,center=0.05,
    cmap="coolwarm",
    square=False
)
ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=90,
    horizontalalignment='right'
)
plt.show()

#miranda
miranda_df=pd.read_csv(miranda_result_tab,index_col=0,sep='\t')
ax = sns.heatmap(
    miranda_df, 
    vmin=0, vmax=1,center=0.05,
    cmap="coolwarm",
    square=False
)
ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=90,
    horizontalalignment='right'
)
plt.show()

#pita
pita_df=pd.read_csv(miranda_result_tab,index_col=0,sep='\t')
ax = sns.heatmap(
    pita_df, 
    vmin=0, vmax=1,center=0.05,
    cmap="coolwarm",
    square=False
)
ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=90,
    horizontalalignment='right'
)
plt.show()


#5)get the top ones and ho in heatmaps

k1 = Counter(g_targetscan_out_info)
high1=k1.most_common(100)
k1_lnc=[];k1_mirna=[]
for i in high1:
	k1_lnc.append(i[0][0])
	k1_mirna.append(i[0][1])
k1_lnc=set(k1_lnc)
k1_mirna=set(k1_mirna)
with open(targetscan_result_out_top100,"w") as targetscan_result_out_top100_write_to:
	targetscan_result_out_top100_write_to.write("{}\n".format("\t".join(k1_lnc)))
	for i in k1_mirna:
		alist=[]
		for j in k1_lnc:
			try:
				alist.append(g_targetscan_out_info[(j,i)])
			except:
				alist.append(0)
		targetscan_result_out_top100_write_to.write("{}\t{}\n".format(i,"\t".join(map(str,alist))))	

targetscan_df_top100=pd.read_csv(targetscan_result_out_top100,index_col=0,sep='\t')
sns.heatmap(
    targetscan_df_top100, 
    vmin=0, vmax=30,center=5,
    cmap="coolwarm",
    square=False,xticklabels=True, yticklabels=True
)
ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=90,
    horizontalalignment='right'
)
plt.show()


k2 = Counter(g_miranda_out_info)
high2=k2.most_common(100)
k2_lnc=[]
k2_mirna=[]
for i in high2:
	k2_lnc.append(i[0][0])
	k2_mirna.append(i[0][1])

k2_lnc=set(k2_lnc)
k2_mirna=set(k2_mirna)
with open(miranda_result_tab_top100,"w") as miranda_result_tab_top100_write_to:
	miranda_result_tab_top100_write_to.write("{}\n".format("\t".join(k2_lnc)))
	for i in k2_mirna:
		alist=[]
		for j in k2_lnc:
			try:
				alist.append(g_miranda_out_info[(j,i)])
			except:
				alist.append(0)
		miranda_result_tab_top100_write_to.write("{}\t{}\n".format(i,"\t".join(map(str,alist))))	

miranda_df_top100=pd.read_csv(miranda_result_tab_top100,index_col=0,sep='\t')

sns.heatmap(
    miranda_df_top100, 
    vmin=0, vmax=500000,center=100,
    cmap="coolwarm",
    square=False,xticklabels=True, yticklabels=True
)
plt.show()

k3 = Counter(g_pita_out_info)
high3=k3.most_common(100)
k3_lnc=[]
k3_mirna=[]
for i in high3:
	k3_lnc.append(i[0][0])
	k3_mirna.append(i[0][1])

k3_lnc=set(k3_lnc)
k3_mirna=set(k3_mirna)

with open(pita_result_tab_top100,"w") as pita_result_tab_top100_write_to:
	pita_result_tab_top100_write_to.write("{}\n".format("\t".join(k3_lnc)))
	for i in k3_mirna:
		alist=[]
		for j in k3_lnc:
			try:
				alist.append(g_pita_out_info[(j,i)])
			except:
				alist.append(0)
		pita_result_tab_top100_write_to.write("{}\t{}\n".format(i,"\t".join(map(str,alist))))	

pita_df_top100=pd.read_csv(pita_result_tab_top100,index_col=0,sep='\t')

sns.heatmap(
    pita_df_top100, 
    vmin=20, vmax=40,center=20,
    cmap="coolwarm",
    square=False,xticklabels=True, yticklabels=True
)
plt.show()
