#!/bin/bash
output_dir=/Users/yiranli/Desktop/polyA/contextmap_result
parent_dir=/Users/yiranli/Desktop/polyA/ContextMap_v2.7.9

#
# Simply modify the following paths to start a ContextMap run with your own data
#
aligner_name="bowtie2"
aligner_bin="/Users/yiranli/Downloads/bowtie2-2.2.9/bowtie2"
indexer_bin="/Users/yiranli/Downloads/bowtie2-2.2.9/bowtie2-build"
contextmap_jar=$parent_dir/ContextMap_v2.7.9.jar
reads="/Users/yiranli/Desktop/0/Crypto_run1_merged.fastq"
##/Users/yiranli/Desktop/0/cpar_24h_hisat2_merged.sort.end1.fastq,/Users/yiranli/Desktop/0/cpar_24h_hisat2_merged.sort.end2.fastq
genomes_dir=/Users/yiranli/Desktop/polyA/reference_sequences
aligner_index=$output_dir/aligner_index/demo_index
gtf="/Users/yiranli/Desktop/cpar.gtf"

#
# building indices with the original indexer
#
mkdir $output_dir
mkdir $output_dir/aligner_index/ >> log.txt 2>&1

if [ "$aligner_name" = bowtie1 -o "$aligner_name" = bowtie2 ]
then ${indexer_bin} -f $genomes_dir/cpar.fa ${aligner_index} >> log.txt 2>&1
else ${indexer_bin} index -p ${aligner_index} $genomes_dir/chr1_ENSG00000198746.fa >> log.txt 2>&1
fi

#
# starting contextmap
#

java -Xms5000m -Xmx5000m -jar "$contextmap_jar" mapper -t 8 -aligner_name $aligner_name -aligner_bin $aligner_bin -indexer_bin $indexer_bin -reads $reads -o $output_dir -indices $aligner_index -genome $genomes_dir -gtf $gtf --polyA >> log.txt 2>&1
echo ""

# for starting a multi threaded run, you have to add the '-t' option. Together with some performance tuning for the virtual machine the call will look like this (using 8 threads):
#
# jar $contextmap_jar mapper -reads $reads -readlen $readlen -o $output_dir -rmap $rmap -bwtbuild $bwtbuild -bwtindex $bwtindex -genome $genome_dir -t 8
