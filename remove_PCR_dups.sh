#!bin/bash

#This script simply takes trimmed Ilumina reads (from trim_reads_paired.sh) and removes PCR duplicates in every pair in a samples list using bbmap and 30Gb of RAM
#to allow for substitutions while filtering duplicates, add subs=[1] after dedupe in each line below
#Made by Brendan J. Pinto on 9/1/2020

#gather samples names from the working directory
ls *_val_*.fq.gz > samples.txt
sed -i "s/_1_val_1.fq.gz//g" samples.txt
sed -i "s/_2_val_2.fq.gz//g" samples.txt
sed -i "s/_R1_val_1.fq.gz//g" samples.txt
sed -i "s/_R2_val_2.fq.gz//g" samples.txt
uniq samples.txt > uniq_samples.txt
rm samples.txt

cat uniq_samples.txt | while read line 
do
echo "dedupe *_[1/2]_val_[1/2].fq.gz..."
clumpify.sh -Xmx30g in=$line\_1_val_1.fq.gz in2=$line\_2_val_2.fq.gz out=$line\_R1.fastq.gz out2=$line\_R2.fastq.gz dedupe 2>/dev/null 
done

cat uniq_samples.txt | while read line 
do
echo "dedupe *_R[1/2]_val_[1/2].fq.gz..."
clumpify.sh -Xmx30g in=$line\_R1_val_1.fq.gz in2=$line\_R2_val_2.fq.gz out=$line\_R1.fastq.gz out2=$line\_R2.fastq.gz dedupe 2>/dev/null 
done


