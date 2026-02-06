#!bin/bash

# This script renames sequences such that a specified integer of the largest X sequences are renamed as "chr_[1-X]", 
# while the rest of the sequences are renamed as "ctg_[1-10000]"
# written 11-2025 -- bjp

# requires samtools and seqkit
# mamba install seqkit samtools -y

# bash chr_namer.sh $1 $2 $3
# $1 == input_genome.fasta[.gz] (original fasta[.gz] file)
# $2 == [number_of_chromosomes] (integer)
# $3 == new_file_name.fasta (renamed fasta file, output will be bgzipped)

# example usage
# bash chr_namer.sh my_genome.fa 23 my_genome_updated.fa

echo "renaming chromosomes..."
seqkit sort -l -r $1 2>/dev/null | awk '/^>/{print ">chr_" ++i; next}{print}' > zzzzzzzzzztmp1 

echo "indexing genome..."
samtools faidx zzzzzzzzzztmp1
head -$2 zzzzzzzzzztmp1.fai | cut -f2 | tail -1 > zzzzzzzzzztmpval
echo "##### (25%)"

echo "splitting reference..."
seqkit seq -m $(< zzzzzzzzzztmpval) zzzzzzzzzztmp1 -o zzzzzzzzzzchrs 2>/dev/null 

seqkit seq -M $(( $(< zzzzzzzzzztmpval) - 1)) $1 -o zzzzzzzzzztmp2   2>/dev/null 
seqkit sort -l -r zzzzzzzzzztmp2 2>/dev/null | awk '/^>/{print ">ctg_" ++i; next}{print}' > zzzzzzzzzzscaf  
echo "########## (50%)"

echo "compiling ouput..."
cat zzzzzzzzzzchrs zzzzzzzzzzscaf > $3
echo "############### (75%)"

echo "bgzipping output..."
rm zzzzzzzzzz*
bgzip -@4 -f $3

echo "printing sequence names for inspection..."
zcat $3\.gz | grep ">"

echo "#################### (100%)"
