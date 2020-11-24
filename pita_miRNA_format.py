mirna_dict={}
with open("/Users/yiranli/Dropbox/master_assembly/miRNA_binding/targetscan_70/miR_Family_Info.txt") as f:
	header=f.readline()
	for line in f:
		features=line.rstrip().split()
		name=features[3]
		seed=features[1]
		species=features[2]
		if species=="9606":
			if name not in mirna_dict.keys():
				mirna_dict[name]=[seed,[species]]
			else:
				if species not in mirna_dict[name][1]:
					mirna_dict[name][1].append(species)
for key,value in mirna_dict.items():
	print('{}\t{}\t{}'.format(key,value[0],','.join(value[1])))