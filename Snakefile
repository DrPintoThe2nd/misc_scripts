import os

configfile: "snps_config.json"

sample = config["samples_ALL"]
raw_reads = config["raw_reads"]
trimmed_reads = config["trimmed_reads"]
cleaned_reads = config["cleaned_reads"]
mapped_reads = config["mapped_reads"]
called_reads = config["called_reads"]

target_genome = config["target_genome"]
target_genome_fai = config["target_genome_fai"]

rule all:
	input:
#fastqc_analysis rule
		expand(os.path.join(raw_reads, "{sample}_1_fastqc.html"), sample=sample),
		expand(os.path.join(raw_reads, "{sample}_2_fastqc.html"), sample=sample),
#trim_galore_pe rule
		expand(os.path.join(trimmed_reads + "{sample}_1_val_1.fq.gz"), sample=sample),
		expand(os.path.join(trimmed_reads + "{sample}_2_val_2.fq.gz"), sample=sample),
#remove_PCR_dups rule
		expand(os.path.join(cleaned_reads + "{sample}_R1.fastq.gz"), sample=sample),
		expand(os.path.join(cleaned_reads + "{sample}_R2.fastq.gz"), sample=sample),
#minimap2_bam rule
		expand(os.path.join(mapped_reads + "{sample}.bam"), sample=sample),
		expand(os.path.join(mapped_reads + "{sample}.bam.bai"), sample=sample),
		expand(os.path.join(mapped_reads + "{sample}.bam.stat"), sample=sample),
#mosdepth rule
		expand(os.path.join(mapped_reads + "{sample}.regions.bed.gz"), sample=sample),
#freebayes_setup rule
		os.path.join(mapped_reads + "bams.txt"),
#freebayes rule
		os.path.join(called_reads + "freebayes_calls.vcf.gz"),

rule fastqc1_analysis:	
	input:
		fqc_in1 = os.path.join(raw_reads + "{sample}_1.fq.gz"),
	output:
		os.path.join(raw_reads + "{sample}_1_fastqc.html"),
	params:
		raw_reads = raw_reads
	shell:
		"""
		fastqc --threads 1 {input.fqc_in1} -o {params.raw_reads};
		"""

rule fastqc2_analysis:	
	input:
		fqc_in2 = os.path.join(raw_reads + "{sample}_2.fq.gz"),
	output:
		os.path.join(raw_reads + "{sample}_2_fastqc.html"),
	params:
		raw_reads = raw_reads
	shell:
		"""
		fastqc --threads 1 {input.fqc_in2} -o {params.raw_reads}
		"""

rule trim_galore_pe:
	input:
		trim = [os.path.join(raw_reads + "{sample}_1.fq.gz"), os.path.join(raw_reads + "{sample}_2.fq.gz")],
	output:
		os.path.join(trimmed_reads + "{sample}_1_val_1.fq.gz"),
		os.path.join(trimmed_reads + "{sample}_1_val_1_fastqc.html"),
		os.path.join(trimmed_reads + "{sample}_2_val_2.fq.gz"),
		os.path.join(trimmed_reads + "{sample}_2_val_2_fastqc.html"),
	params:
		trimmed_reads = trimmed_reads
	shell:
		"""
		trim_galore --paired --cores 4 --fastqc {input.trim} -o {params.trimmed_reads}
		"""

rule remove_PCR_dups:
	input:
		in1 = os.path.join(trimmed_reads + "{sample}_1_val_1.fq.gz"),
		in2 = os.path.join(trimmed_reads + "{sample}_2_val_2.fq.gz"),
	output:
		out1 = os.path.join(cleaned_reads + "{sample}_R1.fastq.gz"),
		out2 = os.path.join(cleaned_reads + "{sample}_R2.fastq.gz"),
	params:
		cleaned_reads = cleaned_reads
	shell:
		"""
		clumpify.sh -Xmx30g dedupe=t tmpdir={params.cleaned_reads} in={input.in1} in2={input.in2} out={output.out1} out2={output.out2} overwrite=t
		"""

rule mosdepth:
	input:
		bam = os.path.join(mapped_reads + "{sample}.bam"),
	output:
		os.path.join(mapped_reads + "{sample}.mosdepth.global.dist.txt"),
		os.path.join(mapped_reads + "{sample}.mosdepth.region.dist.txt"),
		os.path.join(mapped_reads + "{sample}.mosdepth.summary.txt"),
		os.path.join(mapped_reads + "{sample}.regions.bed.gz"),
		os.path.join(mapped_reads + "{sample}.regions.bed.gz.csi"),
	params:
		prefix = os.path.join(mapped_reads + "{sample}"),
	shell:
		"""
		MOSDEPTH_PRECISION=5 mosdepth -x -n -t4 -Q 3 -b 500000 {params.prefix} {input.bam}
		"""

rule minimap2_bam:
	input:
		target = target_genome,
		query1 = os.path.join(cleaned_reads + "{sample}_R1.fastq.gz"),
		query2 = os.path.join(cleaned_reads + "{sample}_R2.fastq.gz"),
	output:
		bam = os.path.join(mapped_reads + "{sample}.bam"),
		bai = os.path.join(mapped_reads + "{sample}.bam.bai"),
		stat = os.path.join(mapped_reads + "{sample}.bam.stat"),
	params:
		param = str("-R '@RG\\tID:{sample}\\tSM:{sample}'")
	shell:
		"""
		minimap2 -ax sr -t6 {params.param} {input.target} {input.query1} {input.query2} | samtools sort -@ 6 - | samtools view -@ 6 -Sb - > {output.bam};
		sambamba index -t4 {output.bam};
		sambamba flagstat -t4 {output.bam} > {output.stat}
		"""

rule freebayes_setup:
	input:
		bams = expand(os.path.join(mapped_reads, "{sample}.bam"), sample=sample),
	output:
		bams = os.path.join(mapped_reads + "bams.txt"),
	params:
		dir = os.path.join(mapped_reads + '*.bam'),
	shell:
		"""
		ls {params.dir} > {output.bams}
		"""

rule freebayes:
	input:
		target = target_genome,
		target_fai = target_genome_fai,
		bam_list = os.path.join(mapped_reads + "bams.txt"),
	output:
		vcf = os.path.join(called_reads + "freebayes_calls.vcf.gz"),
	shell:
		"""
		freebayes-parallel <(fasta_generate_regions.py {input.target_fai} 10000) 8 -f {input.target} -g 200 -L {input.bam_list} | bgzip -c > {output.vcf}
		"""
