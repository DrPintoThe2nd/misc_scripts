#!bin/bash

#This script simply takes paired-end Ilumina reads and runs trim_galore on every pair in a samples list
#Made by Brendan J. Pinto on 9/1/2020

#gather samples names from the working directory
ls *.fq.gz > samples.txt
sed -i "s/_1.fq.gz//g" samples.txt
sed -i "s/_2.fq.gz//g" samples.txt
sed -i "s/_R1.fq.gz//g" samples.txt
sed -i "s/_R2.fq.gz//g" samples.txt
sed -i "s/_1.fastq.gz//g" samples.txt
sed -i "s/_2.fastq.gz//g" samples.txt
sed -i "s/_R1.fastq.gz//g" samples.txt
sed -i "s/_R2.fastq.gz//g" samples.txt
uniq samples.txt > uniq_samples.txt
rm samples.txt

cat uniq_samples.txt | while read line 
do
echo "trimming *_1/2.fq.gz.."
trim_galore --paired $line\_1.fq.gz $line\_2.fq.gz 2>/dev/null 
done

cat uniq_samples.txt | while read line 
do
echo "trimming *_R1/2.fq.gz.."
trim_galore --paired $line\_R1.fq.gz $line\_R2.fq.gz 2>/dev/null 
done 

cat uniq_samples.txt | while read line 
do
echo "trimming *_1/2.fastq.gz.."
trim_galore --paired $line\_1.fastq.gz $line\_2.fastq.gz 2>/dev/null 
done

cat uniq_samples.txt | while read line 
do
echo "trimming *_R1/2.fastq.gz.."
trim_galore --paired $line\_R1.fastq.gz $line\_R2.fastq.gz 2>/dev/null 
done

