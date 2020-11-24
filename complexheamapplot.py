from operator import itemgetter
kmer=[]
kmer_dict={}
x=0
index=0
with open("/Users/yiranli/Dropbox/master_assembly/clusering/4mer_lnc_cds/seekr.csv") as f1:
	header=f1.readline()
	kmer=list(filter(None,header.split(",")))
	x=len(kmer[0])
	for i in kmer:
		index+=1
		A=i.count("A")/x*100
		T=i.count("T")/x*100
		C=i.count("C")/x*100
		G=i.count("G")/x*100
		kmer_dict[index]=[A,T,G,C]

with open("/Users/yiranli/Dropbox/master_assembly/clusering/4mer_lnc_cds/p_order") as f2:
	print("A\tT\tG\tC")
	for k in f2:
		if not k.startswith("x"):
			k=int(k.rstrip())
			print('{}\t{}\t{}\t{}'.format(kmer_dict[k][0],kmer_dict[k][1],kmer_dict[k][2],kmer_dict[k][3]))


