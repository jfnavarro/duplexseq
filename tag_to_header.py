#!/usr/bin/env python
"""
This script extracts UMIs from reads (FASTQ) generated with duplex
sequencing protocol. The UMIs (1 and 2) are extracted from the reads
and append to the read names.

Author: jc.fernandez.navarro@gmail.com
"""

# TODO Add gzip support
# TODO Pass UMI start position as parameter

import sys
import gzip
from argparse import ArgumentParser
from collections import defaultdict

def coroutine(func):
    """
    Coroutine decorator, starts coroutines upon initialization.
    """
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.__next__()
        return cr
    return start

def readfq(fp): # this is a generator function
    """
    Heng Li's fasta/fastq reader function.
    # https://github.com/lh3/readfq/blob/master/readfq.py
    # Unlicensed.
    Parses fastq records from a file using a generator approach.
    :param fp: opened file descriptor
    :returns an iterator over tuples (name,sequence,quality)
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last:break
        #name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

@coroutine
def writefq(fp):  # This is a coroutine
    """
    Fastq writing generator sink.
    Send a (header, sequence, quality) triple to the instance to write it to
    the specified file pointer.
    """
    fq_format = '@{header}\n{sequence}\n+\n{quality}\n'
    try:
        while True:
            record = yield
            read = fq_format.format(header=record[0], sequence=record[1], quality=record[2])
            fp.write(read)
    except GeneratorExit:
        return

def hdr_rename(read_title, read1_tag, read2_tag):
	# This function renames the header with the formatting of
	# *header coordinates,etc*, *tag from read1*, *tag from read2*,
	# *read designation from original header (for paired reads)*
	illumina = read_title.split(" ")[0].split(":")
	if len(illumina) == 7:
		#Illumina CASAVA >=1.8
		#e.g. @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACGcd
		readnum = read_title.split(" ")[1].split(":")[0]
		return "%s|%s%s/%s" % (read_title.split(" ")[0], read1_tag, read2_tag, readnum)
	elif len(illumina) == 5:
		#Illumina CASAVA >=1.4?
		#e.g. @HWUSI-EAS100R:6:73:941:1973#ATCGAT/1
		read_title = read_title.replace(' ', '_')
		return "%s|%s%s/%s" % (read_title.split('/')[0], read1_tag, read2_tag, read_title.split('/')[1])
	else:
		raise ValueError("Unknown read name format: %s" % read_title)

def main():
	parser =  ArgumentParser()
	parser.add_argument('--infile1', dest='infile1', help='Path to FASTQ file for Read 1.', required=True)
	parser.add_argument('--infile2', dest='infile2', help='Path to FASTQ file for Read 2.', required=True)
	parser.add_argument('--taglen', dest='taglen', type=int, default=6,
						help='Length in bases of the duplex tag sequence. [6]')
	parser.add_argument('--spacerlen', dest='spclen', type=int, default=1,
						help='Length in bases of the spacer sequence between duplex tag and the start of target DNA. [1]')
	parser.add_argument('--trimlen', dest='trimlen', type=int, default=0,
						help='Length in bases to trim on both sides of the sequences. [0]')
	o = parser.parse_args()

	readctr = 0
	goodreads = 0
	badtag = 0
	barcode_dict = defaultdict(lambda: 0)

	R1_handle = open(o.infile1, "r")
	R2_handle = open(o.infile2, "r")
	out_R1_handle = open("out_R1.fq", 'w')
	out_R1_writer = writefq(out_R1_handle)
	out_R2_handle = open("out_R2.fq", 'w')
	out_R2_writer = writefq(out_R2_handle)
	for (header1, sequence1, quality1),\
		(header2, sequence2, quality2) in zip(readfq(R1_handle), readfq(R2_handle)):
		readctr += 1
		# UMIs are located at the end of the reads
		tag1, tag2 = sequence1[-o.taglen:], sequence2[-o.taglen:]
		if (tag1.isalpha() and tag1.count('N') == 0) and (tag2.isalpha() and tag2.count('N') == 0):
			out_R1_writer.send((hdr_rename(header1, tag1, tag2),
								sequence1[:-(o.taglen + o.spclen)],
								quality1[:-(o.taglen + o.spclen)]))
			out_R2_writer.send((hdr_rename(header2, tag1, tag2),
								sequence2[:-(o.taglen + o.spclen)],
								quality2[:-(o.taglen + o.spclen)]))
			goodreads += 1
			barcode_dict[tag1 + tag2] += 1
		else:
			badtag += 1
	R1_handle.close()
	R2_handle.close()
	out_R1_handle.close()
	out_R2_handle.close()

	print("Total sequences processed: %s\n" % readctr)
	print("Sequences with passing tags: %s\n" % goodreads)
	print("Bad tags: %s\n" % badtag)

if __name__ == "__main__":
	main()
