import re
file=open("/Users/yiranli/Dropbox/TcI_Brazil_AUGUSTUS_previous_assembly.gff3")
all=""
for line in file:
	if line.startswith("#") and "/" not in line and "-" not in line and any(x.isupper() for x in line) and not any(x.isdigit() for x in line):
		new_line=line.replace("#","").replace(" ","").replace("proteinsequence=[","").strip()
		all+=new_line
all=all.replace("]","\n")
with open("/Users/yiranli/Dropbox/scripts/tmp","w") as f:
	f.write(all)
file=open("/Users/yiranli/Dropbox/scripts/tmp")
with open("/Users/yiranli/Dropbox/scripts/len","w") as write_to:
	for line in file:
		write_to.write("%d\n" % len(line))
