cd 
/Users/yiranli/Desktop/tools/codonW/reader C* -silent 2>error
rm reader.def
#cat error |grep "Warning"|grep -v translatable|grep -E -o "\".*?[[:space:]]"|sed 's/"//g'|sort|uniq >./incomplete.id
#cpar 
cat error |grep "Warning"|grep -v translatable|grep -E -o "[cC]gd\d.*-RA"|sed 's/"//g'|sort|uniq >./incomplete.id 
python /Users/yiranli/Dropbox/scripts/fasta_filter_based_on_name_list.py C* ./incomplete.id
python /Users/yiranli/Dropbox/scripts/combine_cds.py *complete.fa ./cds_combine.fa
cusp ./cds_combine.fa ./cds_combine.cut
head cds_combine.cut

cd /Volumes/Seagate/small_rna/codon_usage/codonw
for i in cpar chom ctyz cubi cmel cand cbai cmur; do cat ${i}/cds_combine.cut|grep -v "#"|awk '{print$1"\t"$4}'|sed '/^$/d' >./${i}.tab; done
join -t $'\t' -a1 -a2 cpar.tab chom.tab |sed '/^$/d' >cpar_chom
join -a1 -a2 -t $'\t' cpar_chom ctyz.tab >cpar_chom_ctyz
join -a1 -a2 -t $'\t'  cpar_chom_ctyz cubi.tab >cpar_chom_ctyz_cubi
join -a1 -a2 -t $'\t'  cpar_chom_ctyz_cubi cmel.tab >cpar_chom_ctyz_cubi_cmel
join -a1 -a2 -t $'\t' cpar_chom_ctyz_cubi_cmel cand.tab >cpar_chom_ctyz_cubi_cmel_cand
join -a1 -a2 -t $'\t' cpar_chom_ctyz_cubi_cmel_cand cbai.tab >cpar_chom_ctyz_cubi_cmel_cand_cbai
join -a1 -a2 -t $'\t' cpar_chom_ctyz_cubi_cmel_cand_cbai cmur.tab >cpar_chom_ctyz_cubi_cmel_cand_cbai_cmur

cd cpar
codonw CryptoDB-42_CparvumIowaII_AnnotatedCDSs.fasta_complete.fa -coa_cu -nomenu -silent 
codonw CryptoDB-42_CparvumIowaII_AnnotatedCDSs.fasta_complete.fa -fop_file fop.coa -cai_file cai.coa -cbi_file cbi.coa  -all_indices -nomenu -silent

cat CryptoDB-42_CparvumIowaII_AnnotatedCDSs.fasta_complete.out|sort -k9,9n|awk '{print$1"\t"$9}'|head -n10|grep -v "titl"|cut -f1 >top10_Nc
cat CryptoDB-42_CparvumIowaII_AnnotatedCDSs.fasta_complete.out|sort -k8,8nr|awk '{print$1"\t"$8}'|head -n10|grep -v "titl"|cut -f1 >top10_fop
cat CryptoDB-42_CparvumIowaII_AnnotatedCDSs.fasta_complete.out|sort -k7,7nr|awk '{print$1"\t"$7}'|head -n10|grep -v "titl"|cut -f1 >top10_cbi
cat CryptoDB-42_CparvumIowaII_AnnotatedCDSs.fasta_complete.out|sort -k6,6nr|awk '{print$1"\t"$7}'|head -n10|grep -v "titl"|cut -f1 >top10_cai