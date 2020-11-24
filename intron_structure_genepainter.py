#!/usr/bin/python3.6
##store 500 intron-containing cpar gene OG##
with open("/Users/yiranli/Dropbox/intron_study/statictics","r") as stats:
	OG_list=[]
	for line in stats:
		try:
			if int(line.strip().split(" ")[9])>0:
				OG_list.append(line.strip().split(" ")[7])
		except ValueError:
			pass	
##build OG-geneID dic#
with open("/Users/yiranli/Dropbox/intron_study/Orthogroups.txt","r") as orthos:
	OG_dict={}
	for line in orthos:
		list=line.replace(':','').strip().split( )
		OG_dict[list[0]]=list[1:]

##extract protein seqs for each OG#
from Bio import SeqIO
proseq_dict=SeqIO.index("/Users/yiranli/Dropbox/intron_study/genepainter/proteins/all10.fa","fasta")
for i in OG_list:
	file='/Users/yiranli/Dropbox/intron_study/genepainter/proteins/OGs/{}.fasta'.format(i)
	with open(file,"a") as output:
		for n in OG_dict[i]:
			try:
				SeqIO.write(proseq_dict[n],output,"fasta")
			except:
				pass
				
				
##build alignment file for each OG protein file#not working, don't know why, run this step on sapelo1
#from Bio.Align.Applications import ClustalwCommandline
#ClustalwCommandline(cmd='clustalw2', infile='OG0000000.fasta', output='fasta')
##PBS -S /bin/bash
##PBS -N j_ClustalW
##PBS -q batch
##PBS -l nodes=1:ppn=1:HIGHMEM
##PBS -l walltime=480:00:00
##PBS -l mem=10gb
   
#for i in *
#do
#/usr/local/apps/clustalw2/latest/bin/clustalw2  -INFILE=$i -ALIGN -TYPE="PROTEIN" -OUTPUT="FASTA" -OUTORDER="ALIGNED"
#done

###extract gff file for each gene and store in separate OG folder
import os
import sys
import re
import gffutils
###read in gff information and stored in db###
db=gffutils.create_db("/Users/yiranli/Dropbox/Functional_Annotation/others/all10.gff",dbfn="all_db",force=True,keep_order=True,merge_strategy='merge', sort_attribute_values=True)
db=gffutils.FeatureDB("all_db",keep_order=True)
OG_list=set(OG_list)
for i in OG_list:
	dir_name='/Users/yiranli/Dropbox/intron_study/genepainter/gene/{}/'.format(i)
	os.makedirs(dir_name)
	for n in OG_dict[i]:
		gene_name=re.search('\w+_\d+',n).group()
		file='{}{}.gff'.format(dir_name,gene_name)
		with open(file,"a") as writeto:
			for m in db.children(gene_name,featuretype=['CDS','mRNA'], order_by='start'):
				writeto.write('{}\n'.format(m))
			
		


		

							
	