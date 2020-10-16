#!/usr/bin/env python
"""
This tool extracts UMIs from pair-end reads (FASTQ) generated with duplex
sequencing protocol. The UMIs (1 and 2) are extracted from the beginning of the reads
(R1 and R2) and are then appended to the read names.

This is a simplified, modified version of the code present in:

https://github.com/Kennedy-Lab-UW/Duplex-Sequencing

Input:
	A pair of FASTQ files in gzip format.

Output:
	The same of FASTQ files in gzip format with the UMIs removed and appended to the header

Author: jc.fernandez.navarro@gmail.com
"""

# TODO allow to pass the start position of the tag

import gzip
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict
from fastqutils import *

def rename_read(read_name, UMI):
	""""
	Appends the UMI given as parameter into the first white space of
	the read name given as parameter.
	Returns the new gene name with the UMI in it.
	"""
	read_id = read_name.split(" ")
	read_id[0] = read_id[0] + "_" + UMI
	return " ".join(read_id)

def main():
	parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('infile1', type=str, help='Path to FASTQ gzipped file for R1')
	parser.add_argument('infile2', type=str, help='Path to FASTQ gzipped file for R2')
	parser.add_argument('--out-prefix', dest='out_prefix', type=str, default="out",
						help='Prefix to prepend to the output files. [out]')
	parser.add_argument('--umi-length', dest='taglen', type=int, default=6,
						help='Length in bases of the UMI sequence. [6]')
	parser.add_argument('--spacer-length', dest='spclen', type=int, default=1,
						help='Length in bases to trim after the UMI [1]')
	parser.add_argument('--filter-same-umi', action="store_true", default=False, dest='sametag_filter',
						help="Discard read pairs whose UMI are identical (A and B)")
	o = parser.parse_args()

	reads_count = 0
	badtags_count = 0
	sametag_count = 0
	barcode_dict = defaultdict(lambda: 0)

	R1_handle = gzip.open(o.infile1, "rt")
	R2_handle = gzip.open(o.infile2, "rt")
	out_R1_handle = gzip.open("{}_R1.fastq.gz".format(o.out_prefix), 'wt')
	out_R1_writer = writefq(out_R1_handle)
	out_R2_handle = gzip.open("{}_R2.fastq.gz".format(o.out_prefix), 'wt')
	out_R2_writer = writefq(out_R2_handle)
	for (header1, seq1, qua1),\
		(header2, seq2, qua2) in zip(readfq(R1_handle), readfq(R2_handle)):
		reads_count += 1
		tag1 = seq1[:o.taglen]
		tag2 = seq2[:o.taglen]
		if o.sametag_filter and tag1 != tag2:
			sametag_count += 1
			continue
		if tag1.isalpha() and tag1.count('N') == 0 and tag2.isalpha() and tag2.count('N') == 0:
			tag = tag1 + tag2
			out_R1_writer.send((rename_read(header1, tag),
								seq1[o.taglen + o.spclen:],
								qua1[o.taglen + o.spclen:]))
			out_R2_writer.send((rename_read(header2, tag),
								seq2[o.taglen + o.spclen:],
								qua2[o.taglen + o.spclen:]))
			barcode_dict[tag] += 1
		else:
			badtags_count += 1
	R1_handle.close()
	R2_handle.close()
	out_R1_handle.close()
	out_R2_handle.close()

	print("Total read pairs processed: {}".format(reads_count))
	print("Unique UMIs found: {}".format(len(barcode_dict.keys())))
	if o.sametag_filter:
		print("Discarded read pairs with same UMI: {}".format(sametag_count))
	print("Discarded read pairs with a bad UMI: {}".format(badtags_count))

if __name__ == "__main__":
	main()
