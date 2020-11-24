####remove host contanimation####
#################################
#1)map against target genome with 4 lane combined data#
bowtie2 -x cpar_DB -1 /Volumes/Seagate/adam_rna-seq/24_and_48hr_HCT8/24hr_2_RNeasy-44994001/24h_S2_R1 -2 /Volumes/Seagate/adam_rna-seq/24_and_48hr_HCT8/24hr_2_RNeasy-44994001/24h_S2_R2  -S 24h_S2_all.sam -p 8
a) bowtie2 mapping against host sequence
Host example: human genome hg19 (download bowtie2 hg19 index)

# 1) create bowtie2 index database (host_DB) from host reference genome 
bowtie2-build host_genome.fna host_DB

# 2) bowtie2 mapping against host sequence database, keep both mapped and unmapped reads (paired-end reads)
bowtie2 -x host_DB -1 SAMPLE_r1.fastq -2 SAMPLE_r2.fastq -S SAMPLE_mapped_and_unmapped.sam

# 3) convert file .sam to .bam
samtools view -bS SAMPLE_mapped_and_unmapped.sam > SAMPLE_mapped_and_unmapped.bam
b) filter required unmapped reads

# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -b -F 12 -F 256 SAMPLE_mapped_and_unmapped.bam > SAMPLE_bothEndsUnmapped.bam
-f 12     Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
-F 256   Do not(-F) extract alignments which are: <not primary alignment>
see meaning of SAM-flags
c)  split paired-end reads into separated fastq files .._r1 .._r2

# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -n SAMPLE_bothEndsUnmapped.bam >SAMPLE_bothEndsUnmapped_sorted 

bedtools bamtofastq -i SAMPLE_bothEndsUnmapped_sorted.bam -fq SAMPLE_host_removed_r1.fastq -fq2 SAMPLE_host_removed_r2.fastq

Result
Two files of paired-end reads, containing non-host sequences
SAMPLE_host_removed_r1.fastq
SAMPLE_host_removed_r2.fastq

###STEP1###########
#####Trimming######

module load java/jdk1.8.0_20 fastqc
time fastqc  /home/yl22225/antisense/v1/Crypto_run1_merged.fastq

module load java/jdk1.8.0_20
java -jar /Users/yiranli/Desktop/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog trimlog /Users/yiranli/Desktop/48h_S3_both_mapped_R1.fastq /Users/yiranli/Desktop/48h_S3_both_mapped_R2.fastq 48h_S3_mapped_R1_trimmed.fastq out1_s3 48h_S3_mapped_R2_trimmed.fastq out2_s3 ILLUMINACLIP:/Users/yiranli/Desktop/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

###STEP2##########
#####Mapping######

hisat2 -p 8 --dta -x /Users/yiranli/Desktop/cpar_tran -1 /Users/yiranli/Desktop/48h_S1_mapped_R1_trimmed.fastq -2 /Users/yiranli/Desktop/48h_S1_mapped_R2_trimmed.fastq --rna-strandness RF --summary-file hisat_summary_48h_S1 -S cpar_48h_S1_hisat2.sam

####bam filter####
##################
#0h
samtools view -b -F 4 -q 10 /Users/yiranli/Desktop/cpar_0h_hisat2.sort.bam -o /Users/yiranli/Desktop/very_cleaned_reads/cpar_0h_hisat2.sort.bam

#24h & 48h

cd /Users/yiranli/Desktop/
for i in cpar_24h_S1_hisat2.sort.bam cpar_24h_S2_hisat2.sort.bam cpar_24h_S3_hisat2.sort.bam cpar_48h_S1_hisat2.sort.bam cpar_48h_S2_hisat2.sort.bam cpar_48h_S3_hisat2.sort.bam 
do 
samtools view -b -f 2 -q 10 $i -o /Users/yiranli/Desktop/very_cleaned_reads/${i}
done

for i in cpar_24h_S1_hisat2.sort.bam cpar_24h_S2_hisat2.sort.bam cpar_24h_S3_hisat2.sort.bam cpar_48h_S1_hisat2.sort.bam cpar_48h_S2_hisat2.sort.bam cpar_48h_S3_hisat2.sort.bam 
do 
samtools sort $i >${i}.sort
done

#####counts#######
##################
python3.6 /Users/yiranli/Desktop/tools/htseq/scripts/htseq-count \
-f bam -t transcript -r pos -s yes -i gene_id \
/Users/yiranli/Desktop/very_cleaned_reads/cpar_0h_S1_hisat2.sort.bam \
/Users/yiranli/Desktop/lncRNA_reverse.gtf \
>/Users/yiranli/Dropbox/antisense/v2/lncRNA_0h_S1_count

for x in 24h 48h
do
for y in S1 S2 S3
do
python3.6 /Users/yiranli/Desktop/tools/htseq/scripts/htseq-count \
-f bam -t transcript -r pos -s reverse -i gene_id /Users/yiranli/Desktop/very_cleaned_reads/cpar_${x}_${y}_hisat2.sort.bam \
/Users/yiranli/Desktop/lncRNA_origin.gtf >/Users/yiranli/Dropbox/antisense/v2/lncRNA_${x}_${y}_count
done
done



###STEP3###########
#####assembly######
stringtie -p 8 -G /home/yl22225/antisense/v1/cpar.gtf -o cpar_0h_transcripts_stringtie.gtf --rf -C -B -l cpar -m 30 -c 5 /home/yl22225/antisense/v1/cpar_0h_hisat2.bam

stringtie /home/yl22225/antisense/v1/cpar_0h_hisat2.bam -p 8 -o cpar_0h_transcripts_stringtie2.gtf --rf -A cpar_0h_gene_abundance.tab -l cpar -m 30 -c 5 


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#after Dec05
cd /Users/yiranli/Desktop/very_cleaned_reads
stringtie -p 8 -G /Users/yiranli/Dropbox/Functional_Annotation/Cpar_BMGF/Cpar_BMGF.gff -o /Users/yiranli/Dropbox/Alternative_splicing/0h_transcript_stringtie.gff -A 0h_transcript_stringtie.abundance --rf -j 5 -c 10 -g 5 -l 0h cpar_0h_S1_hisat2.sort.bam

cd /Users/yiranli/Desktop/very_cleaned_reads
stringtie -p 8 -G /Users/yiranli/Dropbox/Functional_Annotation/Cpar_BMGF/Cpar_BMGF.gff -o /Users/yiranli/Dropbox/Alternative_splicing/24h_transcript_stringtie.gff -A 24h_transcript_stringtie.abundance --rf -j 5 -c 15 -g 5 -l 24h cpar_24h_S1_hisat2.sort.bam cpar_24h_S2_hisat2.sort.bam cpar_24h_S3_hisat2.sort.bam

cd /Users/yiranli/Desktop/very_cleaned_reads
stringtie -p 8 -G /Users/yiranli/Dropbox/Functional_Annotation/Cpar_BMGF/Cpar_BMGF.gff -o /Users/yiranli/Dropbox/Alternative_splicing/48h_transcript_stringtie.gff -A 48h_transcript_stringtie.abundance --rf -j 5 -c 15 -g 5 -l 48h cpar_48h_S1_hisat2.sort.bam cpar_48h_S2_hisat2.sort.bam cpar_48h_S3_hisat2.sort.bam
#-B/-B:  for Ballgown input
#-A: Gene abundances will be reported
#-C: reference file that are fully covered by reads
#-j 2:require at least 2 junction supporting reads
#-c: coverage required for transcript
#g: merge cap

#filter possible runthrough
cd /Users/yiranli/Desktop/very_cleaned_reads
stringtie -p 8 -o \ /Users/yiranli/Dropbox/Alternative_splicing/0h_transcript_stringtie_noguide.gff --rf -l 0h cpar_0h_S1_hisat2.sort.bam

cd /Users/yiranli/Desktop/very_cleaned_reads
stringtie -p 8 -o /Users/yiranli/Dropbox/Alternative_splicing/24h_transcript_stringtie_noguide.gff --rf -l 24h cpar_24h_S1_hisat2.sort.bam cpar_24h_S2_hisat2.sort.bam cpar_24h_S3_hisat2.sort.bam

cd /Users/yiranli/Desktop/very_cleaned_reads
stringtie -p 8 -o \ /Users/yiranli/Dropbox/Alternative_splicing/48h_transcript_stringtie_noguide.gff --rf -l 48h cpar_48h_S1_hisat2.sort.bam cpar_48h_S2_hisat2.sort.bam cpar_48h_S3_hisat2.sort.bam


###separate mRNA from lncRNA by bedtools
for i in 0 24 48
do
bedtools intersect -b /Users/yiranli/Dropbox/Alternative_splicing/cpar_mRNA.gff -a /Users/yiranli/Dropbox/Alternative_splicing/${i}h_transcript_stringtie.gff -wa -f 1.0 -r -e -s>/Users/yiranli/Dropbox/Alternative_splicing/${i}h_coding_transcripts_par.gff 
#some genes only remained partial
##clearly extract from compete gff
#step1 grep id
cat /Users/yiranli/Dropbox/Alternative_splicing/${i}h_coding_transcripts_par.gff|cut -f9|cut -d ";" -f2|sed 's/^ //g'|sed 's/$/;/g'|sort|uniq>/Users/yiranli/Dropbox/Alternative_splicing/${i}h_coding_transcript_id

#step2 get the gff with same transcript id
grep -f /Users/yiranli/Dropbox/Alternative_splicing/${i}h_coding_transcript_id /Users/yiranli/Dropbox/Alternative_splicing/${i}h_transcript_stringtie.gff|sort >/Users/yiranli/Dropbox/Alternative_splicing/${i}h_coding_transcripts.gff

###get the lncRNAs as remaining transcripts after mRNA extraction
grep -v -f /Users/yiranli/Dropbox/Alternative_splicing/${i}h_coding_transcript_id /Users/yiranli/Dropbox/Alternative_splicing/${i}h_transcript_stringtie.gff|sort >/Users/yiranli/Dropbox/Alternative_splicing/${i}h_LNC_transcripts.gff
#remove tRNA from lncRNA
cat /Users/yiranli/Dropbox/Alternative_splicing/${i}h_LNC_transcripts.gff |grep -v 'reference_id "CPIW' >1
mv 1 /Users/yiranli/Dropbox/Alternative_splicing/${i}h_LNC_transcripts.gff
done

####transcript_merge#####
stringtie --merge /Users/yiranli/Dropbox/Alternative_splicing/0h_coding_transcripts.gff /Users/yiranli/Dropbox/Alternative_splicing/24h_coding_transcripts.gff /Users/yiranli/Dropbox/Alternative_splicing/48h_coding_transcripts.gff -G /Users/yiranli/Dropbox/Functional_Annotation/Cpar_BMGF/Cpar_BMGF.gff -i>/Users/yiranli/Dropbox/Alternative_splicing/all_coding_transcripts_intronRetained.gtf

stringtie --merge /Users/yiranli/Dropbox/Alternative_splicing/0h_LNC_transcripts.gff /Users/yiranli/Dropbox/Alternative_splicing/24h_LNC_transcripts.gff /Users/yiranli/Dropbox/Alternative_splicing/48h_LNC_transcripts.gff -i>/Users/yiranli/Dropbox/Alternative_splicing/all_LNC_transcripts_intronRetained.gtf

##filter wired exons with length 1bp##
cd /Users/yiranli/Dropbox/Alternative_splicing/
cat all_coding_transcripts_intronRetained.gtf |awk 'function abs(v) {return v<0? -v:v} {if (abs($4-$5)>2) print}'>1 
mv 1 all_coding_transcripts_intronRetained.gtf
cat all_LNC_transcripts_intronRetained.gtf |awk 'function abs(v) {return v<0? -v:v} {if (abs($4-$5)>2) print}'>1 
mv 1 all_LNC_transcripts_intronRetained.gtf
##then filter transcripts without exons"
grep -f <(cat /Users/yiranli/Dropbox/Alternative_splicing/all_coding_transcripts_intronRetained.gtf |awk '{if ($3=="exon") print}'|cut -f9|cut -d ";" -f1|sort|uniq) all_coding_transcripts_intronRetained.gtf |sort|uniq>1
mv 1 all_coding_transcripts_intronRetained.gtf

##SUPPA
rm -r /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_coding/
mkdir /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_coding/
python3.6 /usr/local/lib/python3.6/site-packages/SUPPA/suppa.py generateEvents -i /Users/yiranli/Dropbox/Alternative_splicing/all_coding_transcripts_intronRetained.gtf -o /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_coding/ -e SE SS MX RI FL
cd /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_coding/
cat *.gtf >all.gtf
cat *.ioe >all.ioe

rm -r /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_LNC/
mkdir /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_LNC/
python3.6 /usr/local/lib/python3.6/site-packages/SUPPA/suppa.py generateEvents -i seqname	gene_id	event_id	alternative_transcripts	total_transcripts
/Users/yiranli/Dropbox/Alternative_splicing/all_LNC_transcripts_intronRetained.gtf -o /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_LNC/ -e SE SS MX RI FL
cd /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_LNC/
cat *.gtf >all.gtf
cat *.ioe >all.ioe
cd /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_coding/
cat *.gtf >all.gtf
cat *.ioe >all.ioe

###prepare file for the report##
##Cpar_gene_withLncRNAOverlapped.gff##
bedtools intersect -a /Users/yiranli/Dropbox/Functional_Annotation/Cpar_BMGF/Cpar_BMGF.gff -b /Users/yiranli/Dropbox/Alternative_splicing/all_LNC_transcripts_intronRetained.gtf|sort|uniq|sort -k4 >/Users/yiranli/Dropbox/BMGF_report_Jan09/new_after_uniqueness_clean/Cpar_gene_withLncRNAOverlapped.gff

##Cpar_LncRNA_candidates.gtf##
cat /Users/yiranli/Dropbox/Alternative_splicing/all_LNC_transcripts_intronRetained.gtf |sed 's/StringTie/./g' >/Users/yiranli/Dropbox/BMGF_report_Jan09/new_after_uniqueness_clean/Cpar_LncRNA_candidates.gtf

##LncRNA_AlternativeSplicing.gtf##
bedtools intersect -b /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_LNC/all.gtf -a /Users/yiranli/Dropbox/Alternative_splicing/all_LNC_transcripts_intronRetained.gtf -wa|sort|uniq|sort -k4>/Users/yiranli/Dropbox/BMGF_report_Jan09/new_after_uniqueness_clean/LncRNA_AlternativeSplicing.gtf

##mRNA_alternativeSplicing.gtf##
bedtools intersect -a /Users/yiranli/Dropbox/Alternative_splicing/all_coding_transcripts_intronRetained.gtf -b /Users/yiranli/Dropbox/Alternative_splicing/SUPPA_coding/all.gtf  -wa|sort -k4|uniq|sed 's/StringTie/./g'>/Users/yiranli/Dropbox/BMGF_report_Jan09/new_after_uniqueness_clean/mRNA_alternativeSplicing.gtf
grep -f <(cat /Users/yiranli/Dropbox/BMGF_report_Jan09/new_after_uniqueness_clean/Cpar_IOWAII_mRNA_AlternativeSplicing_Kissinger.gtf |cut -f9|cut -d ";" -f2) /Users/yiranli/Dropbox/Alternative_splicing/all_coding_transcripts_intronRetained.gtf |sed 's/StringTie/./g' >/Users/yiranli/Dropbox/BMGF_report_Jan09/new_after_uniqueness_clean/Cpar_IOWAII_mRNA_AlternativeSplicing_Kissinger.gtf



###TransDecoder###
/usr/local/apps/transdecoder/4.0.0/util/cufflinks_gtf_genome_to_cdna_fasta.pl all_LNC_transcripts_intronRetained.gtf /home/yl22225/v6/CryptoDB-32_CparvumIowaII_Genome.fasta >lnc_transcripts.fasta

/usr/local/apps/transdecoder/4.0.0/util/cufflinks_gtf_to_alignment_gff3.pl all_LNC_transcripts_intronRetained.gtf >all_LNC_transcripts_intronRetained.gff3

/usr/local/apps/transdecoder/4.0.0/TransDecoder.LongOrfs -S -t lnc_transcripts.fasta 
/usr/local/apps/transdecoder/4.0.0/TransDecoder.Predict -t lnc_transcripts.fasta

/usr/local/apps/transdecoder/4.0.0/util/cdna_alignment_orf_to_genome_orf.pl \
     /home/yl22225/lnc_transcripts.fasta.transdecoder.gff3 \
     all_LNC_transcripts_intronRetained.gff3 \
     lnc_transcripts.fasta > all_LNC_transcripts.fasta.transdecoder.genome.gff3

sed -i 's/*//g' /home/yl22225/lnc_transcripts.fasta.transdecoder.pep
module load interproscan/5.23-62.0 
sh /usr/local/apps/interproscan/5.23-62.0/interproscan.sh -appl Pfam -d ~ -i /home/yl22225/lnc_transcripts.fasta.transdecoder.pep

cat ips_lnc_transcripts.fasta.transdecoder.pep.tsv|cut -f3 -d :|sort|uniq>ips_hit_lnc_name
##remove this ids from all_lnc_transcript, and run to generate related all files##

##statistics##
#lnc overlap with exon
bedtools intersect -a /Users/yiranli/Dropbox/BMGF_report_Jan09/Cpar_LncRNA_candidates.gff -b /Users/yiranli/Desktop/cpar.exon -wa |awk '{if($3=="transcript") print}'|sort|uniq|wc
#protein with alternative splicing#
bedtools intersect -a /Users/yiranli/Desktop/Cpar_V6.gff3 -b /Users/yiranli/Dropbox/BMGF_report_Jan09/new_after_uniqueness_clean/mRNA_alternativeSplicing.gtf -wa|awk '{if ($3=="mRNA") print}'|sort|uniq|wc


##filter lncRNA with length TPM and occurrence
for x in 0h 24h 48h
do
for y in S1 S2 S3
do
stringtie -p 8 -o /Users/yiranli/Dropbox/antisense/v2/mRNA_${x}_${y}.gtf --rf -j 5 -c 30 -g 5 -l ${x}_${y} /Users/yiranli/Desktop/very_cleaned_reads/cpar_${x}_${y}_hisat2.sort.bam
done
done


join \
<(cat /Users/yiranli/Desktop/lncRNA_origin.gtf |cut -f4,5,9|cut -d " " -f1,2|sed 's/ /\t/g'|cut -f1,2,4|sed 's/"//g'|sed 's/;//g'|awk -F " " '{print$3"\t"$1"\t"$2}'|sort|uniq) \ 
<(join \
<(join \
<(join \
<(join \
<(join \
<(join \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_0h_S1_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_24h_S1_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_24h_S2_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_24h_S3_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_48h_S1_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_48h_S2_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_48h_S3_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort))


join \
<(join \
<(join \
<(join \
<(join \
<(join \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_0h_S1_count |grep MSTRG|awk '{print $1"\t"$9}'|sort) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_24h_S1_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_24h_S2_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_24h_S3_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_48h_S1_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_48h_S2_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort)) \
<(cat /Users/yiranli/Dropbox/antisense/v2/lncRNA_48h_S3_abundance |grep MSTRG|awk '{print $1"\t"$9}'|sort))