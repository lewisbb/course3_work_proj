configfile: 'configs/test.yaml'

DATA_DIR = config['data_dir'] 
GENOME_DIR = config['genome_dir']
SAMPLES = config['samples']

rule all:
    input:
        expand('{dd}/mapped/{sample}.Aligned.sortedByCoord.out.bam', dd=DATA_DIR, sample=SAMPLES),
        expand('{dd}/mapped/{sample}.ReadsPerGene.out.tab', dd=DATA_DIR, sample=SAMPLES),
        'genome_indexed_star/',
        f'{DATA_DIR}/raw_multiqc.html',
        f'{DATA_DIR}/mapped_multiqc.html'

rule fastqc_raw:
	input:
		expand('{dd}/raw/{sample}.fastq.gz', dd=DATA_DIR, sample=SAMPLES)
	output:
		directory(f'{DATA_DIR}/raw_fastqc')
	conda:
		'envs/fastqc.yml'
	shell:
		'''
		[[ ! -d {output} ]] && mkdir {output}
		fastqc -q -o {output} {input}
		'''

rule multiqc_raw:
	input:
		f'{DATA_DIR}/raw_fastqc'
	output:
		f'{DATA_DIR}/raw_multiqc.html'
	conda:
		'envs/multiqc.yml'
	shell:
		'multiqc -q --no-data-dir -i "Raw reads quality" -n {output} {input}'

rule trimming_reads:
	input:
		expand('{dd}/raw/{sample}.fastq.gz', dd=DATA_DIR, sample=SAMPLES)
	output:
		expand('{dd}/clean/{sample}.fastq.gz', dd=DATA_DIR, sample=SAMPLES),
	params:
	    outdir = f'{DATA_DIR}/clean/'
	conda:
		'envs/trimming_reads.yml'
	shell:
		'''
		[[ ! -d {params.outdir} ]] && mkdir {params.outdir}
		for sample in {SAMPLES}
		do
			bbduk.sh \
				in={DATA_DIR}/raw/$sample.fastq.gz \
				out={params.outdir}/$sample.fastq.gz \
				ref=adapters
		done
		'''

rule fastqc_clean:
	input:
		expand('{dd}/clean/{sample}.fastq.gz', dd=DATA_DIR, sample=SAMPLES)
	output:
		directory(f'{DATA_DIR}/clean_fastqc')
	conda:
		'envs/fastqc.yml'
	shell:
		'''
		[[ ! -d {output} ]] && mkdir {output}
		fastqc -q -o {output} {input}
		'''

rule STAR_genome_indexing:
	input:
		fasta = f'{GENOME_DIR}/*.fa',
		gtf = f'{GENOME_DIR}/*.gtf'
	output:
		'genome_indexed_star/'
	conda:
		'envs/star.yml'
	shell:
		'''
		star \
			--runMode GenomeGenerate \ 
			--genomeDir {output} \ 
			--genomeFastaFiles {input.fasta} \ 
			--sjdbGTFfile {input.gtf} 
		''' 

rule mapping_reads:
	input:
		fastq = expand('{dd}/clean/{sample}.fastq.gz', dd=DATA_DIR, sample=SAMPLES),
		genome_indexed = 'genome_indexed_star'
	output:
		bam = expand('{dd}/mapped/{sample}.Aligned.sortedByCoord.out.bam', dd=DATA_DIR, sample=SAMPLES),
		counts = expand('{dd}/mapped/{sample}.ReadsPerGene.out.tab', dd=DATA_DIR, sample=SAMPLES),
		outdir = f'{DATA_DIR}/mapped' 
	params:
		in_prefix = f'{DATA_DIR}/clean/',
		gtf = f'{GENOME_DIR}/*.gtf'
	conda:
		'envs/star.yml'
	shell: 	
		'''
		[[ ! -d {params.outdir} ]] && mkdir {params.outdir}
		for sample in {SAMPLES} 
		do
			star \
				--readFilesIn $sample.fastq.gz \
				--readFilesPrefix {params.in_prefix} \
				--readFilesCommand gunzip -c \
				--genomeDir {input.genome_indexed} \
				--sjdbGTFfile {params.gtf} \
				--outFileNamePrefix {params.outdir}/$sample. \
				--outSAMtype BAM SortedByCoordinate \
				--outSAMunmapped Within \
				--outSAMattributes Standard \
				--quantMode GeneCounts
		done
		'''
		
rule multiqc_secondary:
	input:
		f'{DATA_DIR}/clean_fastqc',
		f'{DATA_DIR}/mapped'
	output:
		f'{DATA_DIR}/mapped_multiqc.html'
	conda:
		'envs/multiqc.yml'
	shell:
		'multiqc -q --no-data-dir -i "Mapped reads statistics" -n {output} {input}'