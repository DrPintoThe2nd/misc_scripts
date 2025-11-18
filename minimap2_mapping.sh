#!bin/bash

#script uses bwa-mem2 and samtools to map all paired reads in a directory to the same (indexed) reference genome, then index all bam files for downstream variant calling
#updated for bwa-mem2 by Brendan J. Pinto on 2/5/2021

#bwa-mem2 index genome_file.fasta
#usage "bash bwa_mapping.sh indexed_genome.fasta"

#gather samples names from the working directory
ls *_R1.fastq.gz > samples.txt
sed -i "s/_R1.fastq.gz//g" samples.txt

cat samples.txt | while read line 
do
echo "mapping them reads, king ;)"
minimap2 -ax sr $1 $line\_R1.fastq.gz $line\_R2.fastq.gz -t8 | samtools sort -@ 4 -O bam - -o $line\.bam 2>/dev/null 
done

#cat samples.txt | while read line 
#do
#echo "indexing bam files, king ;)"
#samtools index $line\.bam
#done

