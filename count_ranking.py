import operator
oocyst1_dict={}
oocyst2_dict={}
oocyst3_dict={}
h01_dict={}
h02_dict={}
h03_dict={}
h04_dict={}

length={}
with open("/Users/yiranli/Desktop/paper_prep/cluster/feature_length_lnc") as f:
	for line in f:
		features=line.split()
		lnc_id=features[0]
		lnc_len=features[1]
		length[lnc_id]=int(lnc_len)

with open("/Users/yiranli/Desktop/paper_prep/cluster/lnc27_count.tab") as f1:
	for line in f1:
		if line.startswith("TCONS"):
			features=line.split("\t")
			lnc_id=features[0]
			oocyst1=int(features[11])
			oocyst2=int(features[20])
			oocyst3=int(features[21])
			h01=int(features[13])
			h02=int(features[14])
			h03=int(features[5])
			h04=int(features[12])
			l=length[lnc_id]
			oocyst1_dict[lnc_id]=oocyst1/l;oocyst2_dict[lnc_id]=oocyst2/l;oocyst3_dict[lnc_id]=oocyst3/l
			h01_dict[lnc_id]=h01/l;h02_dict[lnc_id]=h02/l;h03_dict[lnc_id]=h03/l;h04_dict[lnc_id]=h04/l

a=list(enumerate(sorted(oocyst1_dict.items(),key=operator.itemgetter(1),reverse=True)))
ranking_oocyst1={i[1][0]:int(i[0]) for i in a}

b=list(enumerate(sorted(oocyst2_dict.items(),key=operator.itemgetter(1),reverse=True)))
ranking_oocyst2={i[1][0]:int(i[0]) for i in b}

c=list(enumerate(sorted(oocyst3_dict.items(),key=operator.itemgetter(1),reverse=True)))
ranking_oocyst3={i[1][0]:int(i[0]) for i in c}

d=list(enumerate(sorted(h01_dict.items(),key=operator.itemgetter(1),reverse=True)))
ranking_h01={i[1][0]:int(i[0]) for i in d}

e=list(enumerate(sorted(h02_dict.items(),key=operator.itemgetter(1),reverse=True)))
ranking_h02={i[1][0]:int(i[0]) for i in e}

f=list(enumerate(sorted(h03_dict.items(),key=operator.itemgetter(1),reverse=True)))
ranking_h03={i[1][0]:int(i[0]) for i in f}

g=list(enumerate(sorted(h04_dict.items(),key=operator.itemgetter(1),reverse=True)))
ranking_h04={i[1][0]:int(i[0]) for i in g}

ranking={key:ranking_oocyst1[key]+ranking_oocyst2[key]+ranking_oocyst3[key]+ranking_h01[key]+ranking_h02[key]+ranking_h03[key]+ranking_h04[key] for key in ranking_h04.keys()}

result=sorted(ranking.items(),key=operator.itemgetter(1))
for key,value in result:
	print("{}\t{}".format(key,value))


