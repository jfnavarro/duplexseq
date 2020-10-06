# Duplexseq
A pipeline to make consensus sequences from duplex sequencing 

This is a simplified and modified version of the work introduced in 

Duplex Sequencing software Version 3.0 July 11, 2016 Programs by Scott Kennedy(1) (1) 
Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195

https://github.com/Kennedy-Lab-UW/Duplex-Sequencing

# Requirements
* Python 3
* pysam
* biopython
* bwa
* samtools 
* bedtools
* GATK

# Workflow

First, clone the repository and install the pipeline 

´´´
python setup.py install
´´´

Second, extract the tags from the reads and append them to the headers (use --help to see options)

´´´tag_to_header.py R1.fastq.gz R2.fastq.gz´´´

Third, align reads to the reference genome

´´´bwa mem -t THREADS reference.fasta out_R1.fastq.gz out_R2.fastq.gz > aligned.sam´´´

Fourth, sort aligned reads by position

´´´samtools sort --threads THREADS aligned.sam > aligned_sorted.bam´´´

Fifth, create consensus reads (use --help to see options)

´´´create_consensus.py aligned_sorted.bam´´´

Sixth, sort consensus reads by name

´´´samtools sort -n --threads THREADS out.bam > out_sorted_name.bam´´´

Last, convert sorted consensus reads to fastq format (if necessary)

´´´bedtools bamtofastq -i out_sorted_name.bam -fq consensus_reads_R1.fastq -fq2 consensus_reads_R2.fastq´´´


Now you can align with bwa again and compute variants following GATK best practices. 


