#!/bin/bash

#mamba install bedtools seqkit

echo "Finding gaps..."
seqkit locate --bed -m 0 -p "N" $1\.fasta.gz > tmp.bed;

echo "Sorting..."
sort -k1,1 -k2,2n tmp.bed > sort_tmp.bed;

echo "Merging final outputs"
bedtools merge -d 1 -i sort_tmp.bed > $1\.gaps.bed;

echo "Finished! :)"
rm tmp.bed sort_tmp.bed
