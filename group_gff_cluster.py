import os
import re
dict_c_id={}
with open("/Users/yiranli/orange_tmp/cluster_expr_v2.tab") as f:
	for i in f:
		line=i.rstrip().split()
		cluster=line[0]
		ID=line[1]
		dict_c_id[ID]=cluster
with open("/Users/yiranli/Desktop/paper_prep/cluster/final_noncoding_v2.gff","r") as gff:
	with open("/Users/yiranli/Desktop/paper_prep/cluster/v2_C1.gff","w") as f1:
		with open("/Users/yiranli/Desktop/paper_prep/cluster/v2_C2.gff","w") as f2:
			with open("/Users/yiranli/Desktop/paper_prep/cluster/v2_C3.gff","w") as f3:
				with open("/Users/yiranli/Desktop/paper_prep/cluster/v2_C4.gff","w") as f4:
					with open("/Users/yiranli/Desktop/paper_prep/cluster/v2_C5.gff","w") as f5:
						with open("/Users/yiranli/Desktop/paper_prep/cluster/v2_C6.gff","w") as f6:
							with open("/Users/yiranli/Desktop/paper_prep/cluster/v2_C7.gff","w") as f7:
								for line in gff:
									name=re.search(r'TCONS_\d+',line).group()
									try:
										a=dict_c_id[name]
									except:
										pass
									if a=="C1":
										f1.write('{}'.format(line))
									elif a=="C2":
										f2.write('{}'.format(line))
									elif a=="C3":
										f3.write('{}'.format(line))
									elif a=="C4":
										f4.write('{}'.format(line))
									elif a=="C5":
										f5.write('{}'.format(line))
									elif a=="C6":
										f6.write('{}'.format(line))
									elif a=="C7":
										f7.write('{}'.format(line))

									

			
		
