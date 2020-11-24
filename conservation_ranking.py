import sys
import os
import re
os.chdir("/Users/yiranli/Desktop/chromerid/ranking/")
##exonerate2bed_sort
#cat exonerate.out |awk '{if($9=="+") print $6"\t"$7"\t"$8"\t"$2";"$1"\t1\t"$9;else print $6"\t"$8"\t"$7"\t"$2";"$1"\t1\t"$9}'>exonerate.bed
##merge all gff
#cat *.gff >all.gff3
##bedtools
#bedtools intersect -a /Users/yiranli/Desktop/chromerid/ranking/exonerate.bed -b /Users/yiranli/Desktop/chromerid/ranking/all.gff3 -s -wo>./inersect.tab
##make venn
#cat /Users/yiranli/Desktop/chromerid/ranking/inersect.tab |awk -F "\t" '{if($16>200) print}'|grep "CryptoDB-40_Vbrassicafo"|grep -E -o "TCONS_\d+"|sort|uniq >/Users/yiranli/Desktop/chromerid/ranking/venn/express/vbra.txt