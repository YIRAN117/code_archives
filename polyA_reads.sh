#!/bin/bash
#I
#cat /Users/yiranli/Desktop/bwa.result |awk '{if($2==83 && $6~/^[1-9][0-9]M.*S$/ $10~/A{10,}$/) print $0}'
#II
#cat /Users/yiranli/Desktop/bwa.result |awk '{if($2==69 && $10~/T{30,}/) print $0}'
#III
#cat /Users/yiranli/Desktop/bwa.result |awk '{if($2==163 && $6~/^[1-9][0-9]M.*S$/ $10~/A{10,}$/) print $0}'
#I
#cat /Users/yiranli/Desktop/bwa.result |awk '{if($2==99 && $6~/^[1-9]+S.*M/ $10~/^T{10,}/) print $0}'

for i in $(cat /Users/yiranli/Desktop/remove.txt)
do
bioawk -c fastx '{if ($seq~/^TTTTTTTTTTTTTTT/) print ">"$name"\n"$seq}' $i|sed -e 's/^T\{15,\}//g' >2
bowtie2 -x /Users/yiranli/Desktop/cpar_DB -f 2 -S /Users/yiranli/Desktop/polyT.sam -p 8
cat /Users/yiranli/Desktop/polyT.sam |awk '{if($3!="*") print $0}'>>/Users/yiranli/Desktop/polyT_filter.sam
done

cat /Users/yiranli/Desktop/polyT_filter.sam |grep -v "^@"|awk '{if (length($10)>25) print $0}'>/Users/yiranli/Desktop/polyT_filter_25nt.sam
cat /Users/yiranli/Desktop/polyA_filter_25nt.sam |grep -v "@"|awk '{print ">"$1"\n"$10}'|awk '$1~">"{print $0}$1!~">"{tmp="";for(i=1;i<60-length($0)+1;i++){tmp=tmp"N"};print $0""tmp}'>1
samtools view -L intergenic.bed /Users/yiranli/Desktop/polyA_filter_25nt_sorted.bam|awk '{if($6~/^[0-9]+M$/ && $10!~/AAAAAAAAAAAAAAAAAAA/ && $10!~/TTTTTTTTTTTTTTTTTTTTT/) print$1}'>reads.name
samtools view -L intergenic.bed /Users/yiranli/Desktop/polyA_filter_25nt_sorted.bam|awk '{if($6~/^[0-9]+M$/ && $10!~/AAAAAAAAAA/ && $10!~/TTTTTTTTTT/) print$0}'|bioawk -c sam '{if(and($flag,32)) s=reverse($seq);if(and($flag,16)) {s=reverse(revcomp($seq))} print">"$qname"\n"s}'|awk '$1~">"{print $0}$1!~">"{tmp="";for(i=1;i<60-length($0)+1;i++){tmp=tmp"N"};print $0""tmp}'>1
samtools view -L intergenic.bed /Users/yiranli/Desktop/polyA_filter_25nt_sorted.bam|awk '{if($6~/^[0-9]+M$/ && $10!~/AAAAAAAAAAAAA/ && $10!~/TTTTTTTTTTTT/) print$0}'|bioawk -c sam '{if(and($flag,32)) s=reverse($seq);if(and($flag,16)) {s=reverse(revcomp($seq))} print">"$qname"\n"s}'|awk '$1~">"{print $0}$1!~">"{tmp="";for(i=1;i<60-length($0)+1;i++){tmp=tmp"N"};print $0""tmp}'>1