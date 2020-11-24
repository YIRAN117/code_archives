trf_file="/Users/yiranli/Desktop/Y/trf/Y_curated.fa.2.7.7.80.10.50.500.dat"
trf_file_out="/Users/yiranli/Desktop/Y/trf/Y_curated.fa.2.7.7.80.10.50.500.bed"
repeat_masker_file="/Users/yiranli/Desktop/Y/repeatmasker/genome.fasta.out"
repeat_masker_out="/Users/yiranli/Desktop/Y/repeatmasker/genome.fasta.bed"
with open(repeat_masker_file) as f1:
	with open(repeat_masker_out,"w") as f2:
		for line in f1:
			features=line.rstrip().split()
			chro=features[4]
			begin=features[5]
			end=features[6]
			repeat=features[9]
			family=features[10]
			f2.write('{}\t{}\t{}\t{}\t{}\t+\n'.format(chro,begin,end,repeat,family))
a=0
with open(trf_file) as f3:
	with open(trf_file_out,"w") as f4:
		all_text=""
		for line in f3:
			if line.startswith("Sequence") or line[0].isdigit():
				all_text=all_text+str(line)
		all_text=all_text.split("Sequence: ")
		for i in all_text:
			features=list(filter(None,i.split("\n")))
			if len(features) >1:
				chro=features[0]
				for j in features[1:]:
					j_features=list(filter(None,j.split()))
					begin=j_features[0]
					end=j_features[1]
					size=int(j_features[2])
					repeat=j_features[13]
					sequence=j_features[14]
					if size==195:
						a=a+1
						print(">{}\n{}".format(a,repeat))
						f4.write('{}\t{}\t{}\t{}\t{}\t+\n'.format(chro,begin,end,repeat,sequence))

