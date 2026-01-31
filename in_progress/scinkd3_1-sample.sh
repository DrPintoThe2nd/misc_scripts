#!/bin/bash
# define variables
input_ref='GCF_023053635.1_rElgMul1.1.pri_genomic.fna'

male_R1='ELMU_M1_DNAseq_R1.fq.gz'
male_R2='ELMU_M1_DNAseq_R2.fq.gz'
male_meryl='ELMU_M1_DNAseq'

female_R1='ELMU_F1_DNAseq_R1.fq.gz'
female_R2='ELMU_F1_DNAseq_R2.fq.gz'
female_meryl='ELMU_F1_DNAseq'

#source activate scinkd

#count male
time meryl count threads=24 k=28 memory=48 ${male_R1} ${male_R2} output ${male_meryl}\.meryl
#count female
time meryl count threads=24 k=28 memory=48 ${female_R1} ${female_R2} output ${female_meryl}\.meryl

#remove noise
meryl greater-than 1 ${male_meryl}\.meryl   output ${male_meryl}\.1.meryl
meryl greater-than 1 ${female_meryl}\.meryl output ${female_meryl}\.1.meryl

#male-specific
time meryl difference threads=24 ${male_meryl}\.1.meryl ${female_meryl}\.1.meryl output ${male_meryl}\.filt.meryl
#female-specific
time meryl difference threads=24 ${female_meryl}\.1.meryl ${male_meryl}\.1.meryl output ${female_meryl}\.filt.meryl

#meryl_lookup_R1
time meryl-lookup -sequence ${male_R1} -mers ${male_meryl}\.filt.meryl -include -output ${male_meryl}\.filt_R1.fq.gz
time meryl-lookup -sequence ${female_R1} -mers ${female_meryl}\.filt.meryl -include -output ${female_meryl}\.filt_R1.fq.gz
#meryl_lookup_R2
time meryl-lookup -sequence ${male_R2} -mers ${male_meryl}\.filt.meryl -include -output ${male_meryl}\.filt_R2.fq.gz
time meryl-lookup -sequence ${female_R2} -mers ${female_meryl}\.filt.meryl -include -output ${female_meryl}\.filt_R2.fq.gz

#fix pairing issues due to read pair differences
repair.sh -Xmx24g in=${male_meryl}\.filt_R1.fq.gz in2=${male_meryl}\.filt_R2.fq.gz out=${male_meryl}\.filt.paired_R1.fastq.gz out2=${male_meryl}\.filt.paired_R2.fastq.gz overwrite=t
repair.sh -Xmx24g in=${female_meryl}\.filt_R1.fq.gz in2=${female_meryl}\.filt_R2.fq.gz out=${female_meryl}\.filt.paired_R1.fastq.gz out2=${female_meryl}\.filt.paired_R2.fastq.gz outs=tmp overwrite=t

#alignment
minimap2 -t12 -ax sr ${input_ref} ${male_meryl}\.filt.paired_R1.fastq.gz ${male_meryl}\.filt.paired_R2.fastq.gz | samtools sort --write-index -@4 -O bam -o ${male_meryl}\.bam##idx##${male_meryl}\.bam.bai
minimap2 -t12 -ax sr ${input_ref} ${female_meryl}\.filt.paired_R1.fastq.gz ${female_meryl}\.filt.paired_R2.fastq.gz | samtools sort --write-index -@4 -O bam -o ${female_meryl}\.bam##idx##${female_meryl}\.bam.bai
#read depth
mosdepth -Q 10 -x -n -t12 --by 500000 ${male_meryl} ${male_meryl}\.bam
mosdepth -Q 10 -x -n -t12 --by 500000 ${female_meryl} ${female_meryl}\.bam
