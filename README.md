# Duplexseq
A pipeline to make consensus reads from duplex sequencing with UMIs

This is a simplified and modified version of the workflow introduced in:

Duplex Sequencing software Version 3.0 July 11, 2016 Programs by Scott Kennedy(1) (1) 
Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195

In summary, UMIs are extracted from the reads and appended to the headers. 
The reads are then clustered by tag-position and the consensus reads are created
by collapsing these, finally duplex are made using the consensus reads. 

# Notes
- UMIs must be located at the beginning of the reads (0 position)
- Different filters can be applied when creating the consensus reads (check --help)

# Requirements
* Python 3
* pysam
* cutadapt
* bwa
* samtools 
* GATK

# Workflow

Clone the repository and install the pipeline (only performed one time)

```python setup.py install```

First, remove adapters if any

```cutadapt -a ADAPT1 -A ADAPT2 --action mask -j 0 -o R1.fastq.gz -p R2.fastq.gz in1.fastq.gz in2.fastq.gz```

Second, extract the tags from the reads and append them to the headers (use --help to see options)

```extract_umis.py R1.fastq.gz R2.fastq.gz```

Third, align reads to the reference genome

```bwa mem -t THREADS reference.fasta out_R1.fastq.gz out_R2.fastq.gz > aligned.sam```

Fourth, sort aligned reads by position

```samtools sort --threads THREADS aligned.sam > aligned_sorted.bam```

Fifth, create consensus reads (use --help to see options)

```create_consensus.py [options] aligned_sorted.bam```

This tool will generate three pairs of FASTQ files:

- consensus_R1/2 = all the reads (unpaired and paired) where 
a consensus was made (number of reads used in the consensus 
is included in the reads names).

- duplex_R1/2 = all the consensus reads where a duplex was found
and made in both pairs

- duplex_unpaired_R1/2 = all the consensus reads where a duplex
was found and made in either one of the pairs.

Now you can align with bwa again and compute variants following GATK best practices 
(remember to not perform the MarkDuplicates step). 

# License
MIT License, see LICENSE file.

# Authors
Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>



