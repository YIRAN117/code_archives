blastn -db ../genomes/all_genome.fa -query ../../homology/rfam_scan_seq_unmapped.fa -outfmt 5 -evalue 1e-5 -task blastn  >rfam_unmapped.xml
python ~/bin/xml2bed.py -c ext -o ./rfam_unmapped_crypto.tab ./rfam_unmapped.xml

#cat ../../mapping/genome.fasta |sed 's/Contig/cpar/g'|sed 's/.1$/1 | organism=Cryptosporidium_parvum/g' >./cpar.fa
#blastn -db ../../mapping/genome.fasta -query ./unmapped.fa -outfmt 5 -evalue 1e-10 >cpar.xml
#python ~/bin/xml2bed.py -c ext -o cpar.tab cpar.xml
#cat cpar.tab |sed 's/^cpar/Contig/g' >1
#mv 1 cpar.tab
#cat cpar.tab >>rfam_unmapped_crypto.tab


cat rfam_unmapped_crypto.tab |sort -k1,1 -k2,2 -k12,12n >1
mv 1 rfam_unmapped_crypto.tab
cat rfam_unmapped_crypto.tab |cut -f1|sort|uniq >seq.list
while read i; do cat ./rfam_unmapped_crypto.tab|grep "$i" >./result_split_by_contig/${i}.tab; done < seq.list
cd result_split_by_contig/
for i in Con*; do cat $i|awk -F "\t" '{print ">"$25"\n"$22}'|sed 's/| organism=//g'|sed 's/|.*//g'|sed 's/-//g'> ../align/${i%tab}fa; done

for i in *; do mafft --maxiterate 1000 --globalpair ./"${i}" >${i%.fa}_align.fa; done