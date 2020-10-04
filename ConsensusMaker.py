#!/usr/bin/env python
"""
This tools creates consensus reads using aligned pired-end reads with duplex tags
in their headers.

This is a simplified, modifed version of the code present in:

Inputs: 
	A position-sorted paired-end BAM file containing reads with a duplex tag in the header.

Outputs:
	1: A BAM file containing SSCSs
	2: A BAM file containing discarded reads

	Note that quality scores in outputs 1, 2, and 3 are just space fillers and do not signify anything about the
	quality of the sequence.

Author: jc.fernandez.navarro@gmail.com
"""
import pysam
from collections import defaultdict, Counter
from argparse import ArgumentParser

# TODO add option to clusters tags allowing for mismatches
# TODO add option to cluster tags with a window size in position
# TODO add option to clusters taking into consideration the cigar string (mapping end coordinate)
# TODO add option to allow for reads with different lengths

def consensus_maker(grouped_reads_list, cut_off):
    # Reads must have the same length
    # Create a consensus read using the most_common letter from each position
    read_length = len(grouped_reads_list[0])
    num_reads = float(len(grouped_reads_list))
    counters = [Counter([x[i] for x in grouped_reads_list]) for i in xrange(read_length)]
    consensus_read = ''.join(
        [c.most_common()[0][0] if (c.most_common()[0][1] / num_reads) > cut_off else 'N' for c in counters])
    return consensus_read

def main():
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="output BAM file", required=True)
    parser.add_argument('--min_reads', type=int, default=3, dest='min_reads',
                        help="Minimum number of reads allowed to comprise a consensus. [3]")
    parser.add_argument('--max_reads', type=int, default=1000, dest='max_reads',
                        help="Maximum number of reads allowed to comprise a consensus. [1000]")
    parser.add_argument('--cut_off', type=float, default=.7, dest='cut_off',
                        help="Percentage of nucleotides at a given position in a read that must be identical in order "
                             "for a consensus to be called at that position. [0.7]")
    parser.add_argument('--filter-soft-clip', type=str, action="store_true", default=False, dest='softclip_filter',
                        help="Discard reads that are soft-clipped in their alignment.")
    parser.add_argument('--filter-overlap', type=str, action="store_true", default=False, dest='overlap_filter',
                        help="Discard reads overlapping with their pair-end mate.")
    parser.add_argument('--filter-pair', type=str, action="store_true", default=False, dest='pair_filter',
                        help="Discard reads that do not have a proper pair-end mate.")
    parser.add_argument('--filter-singleton', type=str, action="store_true", default=False, dest='single_filter',
                        help="Discard reads that are singletons (only one pair is aligned).")
    parser.add_argument('--max_N', type=float, default=1, dest='Ncut_off',
                        help="Maximum fraction of Ns allowed in a consensus [1.0]")
    o = parser.parse_args()

    # Variables
    reads_count = 0
    reads_count_filtered = 0
    discarded_count = 0
    consensus_count = 0
    read_dict = defaultdict(lambda: defaultdict(list))
    consensus_dict = {}
    min_reads = o.min_reads
    max_reads = o.max_reads
    homo_count = o.rep_filt
    cut_off = o.cut_off
    Ncut_off = o.Ncut_off
    do_overlap_filter = o.softclip_filter
    do_softclip_filter = o.overlap_filter
    do_pair_filter = o.pair_filter
    do_single_filter = o.single_filter
    do_N_filter = Ncut_off < 1.0 and Ncut_off > 0.0

    # Start going through the input BAM file, one position at a time to group reads by start_pos + tag
    in_bam_file = pysam.Samfile(o.infile, "rb")
    discarded_bam_file = pysam.Samfile('discarded.bam', "wb", template=in_bam_file)
    out_bam_file = pysam.Samfile(o.outfile, "wb", template=in_bam_file)
    for record in in_bam_file.fetch(until_eof=True):
        #  Only consider mapped reads
        if record.is_unmapped:
            continue

        reads_count += 1

        #  Get tag
        tag = record.query_name.split('|')[1].split('/')[0]
        tag += ":1" if record.is_read1 else ":2"

        # Mapping filters
        singleton = do_single_filter and record.mate_is_unmapped
        proper_pair = do_pair_filter and record.is_proper_pair

        # Overlap filter
        overlap = False
        if do_overlap_filter and proper_pair and not singleton:
            overlap = (record.reference_start <= record.next_reference_start < (record.reference_start + record.query_length)) or \
                      (record.next_reference_start <= record.reference_start < (record.next_reference_start + record.query_length))

        # soft_clip filter
        soft_clip = False
        if do_softclip_filter:
            for cigar_tuple in record.cigartuples if record.cigartuples else []:
                if cigar_tuple[0] == 4:
                    soft_clip = True
                    break

        # Compute homopolymers count
        homo_A = tag.count('A') >= homo_count
        homo_C = tag.count('C') >= homo_count
        homo_G = tag.count('G') >= homo_count
        homo_T = tag.count('T') >= homo_count

        # Check if the given read is good
        if any([overlap, soft_clip, homo_A, homo_C, homo_G, homo_T, singleton, proper_pair]):
            discarded_count += 1
            discarded_bam_file.write(record)
        else:
            reads_count_filtered += 1
            read_dict[record.reference_start][tag].append(record)

    in_bam_file.close()
    discarded_bam_file.close()

    # Iterate grouped reads to build consensus reads
    for pos, tag_records in read_dict.items():
        for tag, records in tag_records.items():
            sequences = [x.query_sequence for x in records]
            num_sequences = len(sequences)
            if min_reads <= num_sequences <= max_reads:
                #  All sequences same length filter
                first_length = len(sequences[0])
                if not all(len(x) == first_length for x in sequences):
                    continue
                consensus = consensus_maker(sequences, cut_off)
                # Filter out consensuses with too many Ns in them
                if do_N_filter and (consensus.count("N") / float(len(consensus)) >= Ncut_off):
                    continue
                consensus_count += 1
                # Create SAM record (assuming all records in the same position share common attributes)
                record = records[0]
                a = pysam.AlignedRead()
                a.query_name = tag + ":" + str(num_sequences)
                a.flag = record.flag
                a.query_sequence = consensus
                a.reference_id = record.reference_id
                a.reference_start = record.reference_start
                a.mapping_quality = 255
                a.cigartuples = record.cigartuples
                a.next_reference_id = record.next_reference_id
                a.next_reference_start = record.next_reference_start
                a.template_length = record.template_length
                a.qual = 'J' * record.query_length
                # Write SAM record
                altTag = tag.replace(("1" if "1" in tag else "2"), ("2" if "1" in tag else "1"))
                if altTag in consensus_dict:
                    if a.is_read1:
                        out_bam_file.write(a)
                        out_bam_file.write(consensus_dict.pop(altTag))
                    else:
                        out_bam_file.write(consensus_dict.pop(altTag))
                        out_bam_file.write(a)
                else:
                    consensus_dict[tag] = a

    out_bam_file.close()

    # Write summary statistics
    print("Summary Statistics:")
    print("Reads processed: {}".format(reads_count))
    print("Reads discarded: {}".format(discarded_count))
    print("Reads passing filters: {}".format(reads_count_filtered))
    print("Consensus Made: {}".format(consensus_count))

if __name__ == "__main__":
    main()
