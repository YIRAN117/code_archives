import os
import sys
import re
import gffutils
import pybedtools
from pybedtools import BedTool
os.chdir("/Users/yiranli/Desktop/intron")
os.getcwd()
db = gffutils.create_db("/Users/yiranli/Desktop/paper_prep/lncRNA_add_intron.gff",dbfn="cpar_db",force=True,keep_order=True,merge_strategy='merge', sort_attribute_values=True)
db = gffutils.FeatureDB("cpar_db",keep_order=True)
#db.update(db.create_introns(exon_featuretype='exon', grandparent_featuretype='gene', new_featuretype='intron', merge_attributes=True))
#gene=db["cgd8_3173-RA-p1.1"]
#[f.start for f in db.children("cgd8_3173-RA-p1.1",featuretype='CDS')]
with open("/Users/yiranli/Desktop/intron/intron_donor_lnc.bed","w") as f1:
	with open("/Users/yiranli/Desktop/intron/intron_acceptor_lnc.bed","w") as f2:
		for i in db.features_of_type('intron'):
			if i.strand=="+":
				donor_a=int(i.start)-11
				donor_b=int(i.start)+9
				acceptor_a=int(i.end)-21
				acceptor_b=int(i.end)+9
				f1.write('{}\t{}\t{}\t{}\t1\t{}\n'.format(i.seqid,donor_a,donor_b,i.id,i.strand))
				f2.write('{}\t{}\t{}\t{}\t1\t{}\n'.format(i.seqid,acceptor_a,acceptor_b,i.id,i.strand))
			else:
				donor_a=int(i.end)+10
				donor_b=int(i.end)-10
				acceptor_a=int(i.start)+20
				acceptor_b=int(i.start)-10
				f1.write('{}\t{}\t{}\t{}\t1\t{}\n'.format(i.seqid,donor_b,donor_a,i.id,i.strand))
				f2.write('{}\t{}\t{}\t{}\t1\t{}\n'.format(i.seqid,acceptor_b,acceptor_a,i.id,i.strand))
fasta = pybedtools.example_filename('/Users/yiranli/Desktop/CryptoDB-32_CparvumIowaII_Genome.fasta')
donor = pybedtools.BedTool("/Users/yiranli/Desktop/intron/intron_donor_lnc.bed")
donor = donor.sequence(fi=fasta,s=True,name=True).save_seqs("/Users/yiranli/Desktop/intron/intron_donor_lnc.fa") 
acceptor = pybedtools.BedTool("/Users/yiranli/Desktop/intron/intron_acceptor_lnc.bed")
acceptor = acceptor.sequence(fi=fasta,s=True,name=True).save_seqs("/Users/yiranli/Desktop/intron/intron_acceptor_lnc.fa")







