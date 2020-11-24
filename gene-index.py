###make index for genes based on position####
with open("/Users/yiranli/Dropbox/gene.gff","r") as genegff:
	index=0
	dict_gene={}
	for line in genegff:
			list=line.strip().split("\t")
			chr=list[0]
			start=list[3]
			end=list[4]
			strand=list[6]
			index=index+1
			gene_id=list[8].strip().split(";")[0].replace("gene_id=","")
			dict_gene[index]=['{} {} {} {} {}'.format(gene_id,start,end,strand,chr)]
	print(dict_gene[74])
			
####find corresponding sense gene based on antisense start site####
with open("/Users/yiranli/Dropbox/0h_antisense.final.gff") as TPM0:
	a=0
	b=0
	for line in TPM0:
		list=line.strip().split("\t")
		chr=list[1]
		strand=list[4]
		antisense_id=list[0]
		b=b+1
		print(b)
		for keys,values in dict_gene.items():
			if keys==4018:
				continue
			A=str(values).strip().split(" ")[1]#start#
			B=str(values).strip().split(" ")[2]#end#
			C=str(values).strip().split(" ")[4].replace("']","")#chr#
			###find related genes###
			if strand=="+":
				TS=list[2]
				TE=list[3]
				if (int(TS)+100) < int(B) and (int(TS)+100) >int(A) and C==chr:
					sense_id=str(dict_gene[keys]).strip().split(" ")[0].replace("['","")
					sense_strand=str(dict_gene[keys]).strip().split(" ")[3]
					sense_start=str(dict_gene[keys]).strip().split(" ")[1]
					sense_end=str(dict_gene[keys]).strip().split(" ")[2]
					upstream_id=str(dict_gene[int(keys)-1]).strip().split(" ")[0].replace("['","")
					upstream_strand=str(dict_gene[int(keys)-1]).strip().split(" ")[3]
					upstream_start=str(dict_gene[int(keys)-1]).strip().split(" ")[1]
					upstream_end=str(dict_gene[int(keys)-1]).strip().split(" ")[2]
					downstream_id=str(dict_gene[int(keys)+1]).strip().split(" ")[0].replace("['","")
					downstream_strand=str(dict_gene[int(keys)+1]).strip().split(" ")[3]
					downstream_start=str(dict_gene[int(keys)+1]).strip().split(" ")[1]
					downstream_end=str(dict_gene[int(keys)+1]).strip().split(" ")[2]
					a=a+1
					with open("/Users/yiranli/Desktop/0h_related_gene ","a") as write_to:
						write_to.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(a,antisense_id,strand,TS,TE,sense_id,sense_strand,sense_start,sense_end,upstream_id,upstream_strand,upstream_start,upstream_end,downstream_id,downstream_strand,downstream_start,downstream_end))
			else:
				TS=list[3]
				TE=list[2]
				if (int(TS)-100) < int(B) and (int(TS)-100) >int(A) and C==chr:
					sense_id=str(dict_gene[keys]).strip().split(" ")[0].replace("['","")
					sense_strand=str(dict_gene[keys]).strip().split(" ")[3]
					sense_start=str(dict_gene[keys]).strip().split(" ")[1]
					sense_end=str(dict_gene[keys]).strip().split(" ")[2]
					upstream_id=str(dict_gene[int(keys)+1]).strip().split(" ")[0].replace("['","")
					upstream_strand=str(dict_gene[int(keys)+1]).strip().split(" ")[3]
					upstream_start=str(dict_gene[int(keys)+1]).strip().split(" ")[1]
					upstream_end=str(dict_gene[int(keys)+1]).strip().split(" ")[2]
					downstream_id=str(dict_gene[int(keys)-1]).strip().split(" ")[0].replace("['","")
					downstream_strand=str(dict_gene[int(keys)-1]).strip().split(" ")[3]
					downstream_start=str(dict_gene[int(keys)-1]).strip().split(" ")[1]
					downstream_end=str(dict_gene[int(keys)-1]).strip().split(" ")[2]
					a=a+1
					with open("/Users/yiranli/Desktop/0h_related_gene ","a") as write_to:
						write_to.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(a,antisense_id,strand,TS,TE,sense_id,sense_strand,sense_start,sense_end,upstream_id,upstream_strand,upstream_start,upstream_end,downstream_id,downstream_strand,downstream_start,downstream_end))
	
		