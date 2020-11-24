from typing import Dict, List

from Bio import SeqIO
import re

OG_dict: Dict[str, List[str]] = {}
with open("/Users/yiranli/Downloads/Othologs_Crypto/Orthogroups.txt") as f1:
    for line in f1:
        entry = line.rstrip().split()
        OG = entry[0].replace(":", "")
        gene_ids = entry[1:]
        gene_ids: List[str] = [i.replace("-p1", "") for i in gene_ids]
        OG_dict[OG] = gene_ids

cpar_name_dict={}
with open("/Users/yiranli/Downloads/Othologs_Crypto/CryptoDB-46_CparvumIowaII_AnnotatedTranscripts.fasta") as f2:
    for line in f2:
        if line.startswith(">"):
            try:
                cpar_gene_id=re.search("[Cc]gd\d_\d+-RA",line).group()
                cpar_gene_name=re.search("gene_product=.*?\|",line).group().replace("gene_product=","").replace(" |","")
            except:
                print("error:",line)
            cpar_name_dict[cpar_gene_id]=cpar_gene_name



seq_dict = SeqIO.to_dict(SeqIO.parse("/Users/yiranli/Downloads/Othologs_Crypto/all.fa", "fasta"))
print("OG\tname\tcfel\tcpar\tchom\tcmel\tcmur\tcubi\tcand")
for key, value in OG_dict.items():
    cfel = []
    cpar = []
    cand = []
    cmel = []
    cmur = []
    cubi = []
    chom = []
    for i in value:
        if "C_baileyi" not in i and "C_chipmunk" not in i:
            try:
                i_seq=seq_dict[i].seq
            except KeyError:
                print("key error:",i)
            i_seq_GC = str((i_seq.count("G") + i_seq.count("C")) / len(i_seq) * 100)
            if "CMU_" in i:
                cmur.append(i_seq_GC)
            elif "Cfel" in i:
                cfel.append(i_seq_GC)
            elif "cgd" in i:
                try:
                    cpar_name=cpar_name_dict[i]
                except:
                    print("error: name not found",i)
                cpar.append(i_seq_GC)
            elif "cand" in i:
                cand.append(i_seq_GC)
            elif "CmeUK" in i:
                cmel.append(i_seq_GC)
            elif "cubi" in i:
                cubi.append(i_seq_GC)
            elif "Chro." in i:
                chom.append(i_seq_GC)
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(key,cpar_name,",".join(cfel),",".join(cpar),",".join(chom),",".join(cmel),",".join(cmur),",".join(cubi),",".join(cand)))

