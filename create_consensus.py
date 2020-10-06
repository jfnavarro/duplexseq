#!/usr/bin/env python
"""
This tool creates consensus reads using aligned paired-end reads with duplex tags
in their headers.

This is a simplified, modified version of the code present in:

https://github.com/Kennedy-Lab-UW/Duplex-Sequencing

Input:
	A position-sorted paired-end BAM file containing reads with a duplex tag in the header.

Output:
	1: A BAM file containing the duplexs
	2: A BAM file containing discarded reads

Author: jc.fernandez.navarro@gmail.com
"""
import pysam
from Bio.Seq import Seq
from collections import defaultdict, Counter
from argparse import ArgumentParser

# TODO add option to clusters tags allowing for mismatches
# TODO add option to cluster tags with a window size in position

def consensus_maker(grouped_reads_list, cut_off):
    # Reads must have the same length
    # Create a consensus read using the most_common letter from each position
    read_length = len(grouped_reads_list[0])
    num_reads = float(len(grouped_reads_list))
    counters = [Counter([x[i] for x in grouped_reads_list]) for i in range(read_length)]
    consensus_read = ''.join(
        [c.most_common()[0][0] if (c.most_common()[0][1] / num_reads) > cut_off else 'N' for c in counters])
    return consensus_read

def consensus_quality(qual_list):
    # Qualities must have the same length
    # Create a consensus quality using the average of the qualities
    qual_avg = [int(sum(qual_score) / len(qual_list)) for qual_score in zip(*qual_list)]
    return [x if x < 41 else 41 for x in qual_avg]

def main():
    parser = ArgumentParser()
    parser.add_argument('infile', type=str, help="input BAM file")
    parser.add_argument("--outfile", action="store", dest="outfile", help="output BAM file [out.bam]",
                        default="out.bam")
    parser.add_argument('--min-reads', type=int, default=3, dest='min_reads',
                        help="Minimum number of reads allowed to comprise a consensus. [3]")
    parser.add_argument('--max-reads', type=int, default=1000, dest='max_reads',
                        help="Maximum number of reads allowed to comprise a consensus. [1000]")
    parser.add_argument('--homo-count', type=int, default=6, dest='rep_filt',
                        help="Maximum number of homopolymer supported in the tag [A, C, G, T]. [6]")
    parser.add_argument('--cutoff', type=float, default=.7, dest='cut_off',
                        help="Percentage of nucleotides at a given position in a read that must be identical in order "
                             "for a consensus to be called at that position. [0.7]")
    parser.add_argument('--filter-soft-clip', action="store_true", default=False, dest='softclip_filter',
                        help="Discard reads that are soft-clipped in their alignment.")
    parser.add_argument('--filter-overlap', action="store_true", default=False, dest='overlap_filter',
                        help="Discard reads overlapping with their pair-end mate.")
    parser.add_argument('--filter-pair', action="store_true", default=False, dest='pair_filter',
                        help="Discard reads that do not have a proper pair-end mate.")
    parser.add_argument('--filter-singleton', action="store_true", default=False, dest='single_filter',
                        help="Discard reads that are singletons (only one pair is aligned).")
    parser.add_argument('--maxN', type=float, default=1, dest='Ncut_off',
                        help="Maximum fraction of Ns allowed in a consensus [1.0]")
    o = parser.parse_args()

    # Variables
    reads_count = 0
    reads_count_filtered = 0
    discarded_count = 0
    consensus_count = 0
    duplex_count = 0
    read_dict = defaultdict(lambda: defaultdict(list))
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
    print("Grouping reads by position-tag...")
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
        proper_pair = do_pair_filter and not record.is_proper_pair

        # Overlap filter
        overlap = False
        if do_overlap_filter and proper_pair and not singleton:
            overlap = (record.reference_start <= record.next_reference_start < (
                    record.reference_start + record.query_length)) or \
                      (record.next_reference_start <= record.reference_start < (
                              record.next_reference_start + record.query_length))

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
    print("Creating consensus duplexs...")
    for pos, tag_records in read_dict.items():
        consensus_dict = {}
        for tag, records in tag_records.items():
            sequences = [x.query_sequence for x in records]
            qualities = [x.query_qualities for x in records]
            num_sequences = len(sequences)
            if num_sequences >= min_reads and num_sequences <= max_reads:
                #  All sequences same length filter
                first_length = len(sequences[0])
                if not all(len(x) == first_length for x in sequences):
                    continue
                consensus = consensus_maker(sequences, cut_off)
                consensus_qual = consensus_quality(qualities)
                # Filter out consensuses with too many Ns in them
                if do_N_filter and (consensus.count("N") / float(len(consensus)) >= Ncut_off):
                    continue
                consensus_count += 1
                # Create SAM record (assuming all records in the same position share common attributes)
                record = records[0]
                a = pysam.AlignedRead()
                a.query_name = tag
                a.flag = record.flag
                a.query_sequence = consensus
                a.reference_id = record.reference_id
                a.reference_start = record.reference_start
                a.mapping_quality = 255
                a.cigartuples = record.cigartuples
                a.next_reference_id = record.next_reference_id
                a.next_reference_start = record.next_reference_start
                a.template_length = record.template_length
                a.query_qualities = consensus_qual
                #  Store records (there must be one record per tag)
                consensus_dict[tag] = a
        # Iterate consensus dict to find duplexs otherwise write simplex
        processed = set()
        for tag, record in consensus_dict.items():
            #  Check if pair has been processed already
            if tag in processed:
                continue
            clean_tag = tag.split(":")[0]
            index = "2" if record.is_read1 else "1"
            switch_tag = "{}{}:{}".format(clean_tag[int(len(clean_tag) / 2):],
                                          clean_tag[:int(len(clean_tag) / 2)],
                                          index)
            try:
                simplex1 = record.query_sequence
                simplex2 = consensus_dict[switch_tag].query_sequence
                #  Add to processed
                processed.add(tag)
                processed.add(switch_tag)
                if len(simplex1) != len(simplex2):
                    print("Found duplex pair with different lengths\n{}\n{}".format(simplex1, simplex2))
                    continue
                duplex = consensus_maker([simplex1, simplex2], 1.0)
                # Filter out duplex with too many Ns in them
                if do_N_filter and (duplex.count("N") / float(len(duplex)) >= Ncut_off):
                    continue
                duplex_count += 1
            except KeyError:
                duplex = record.query_sequence
            # Create SAM record
            a = pysam.AlignedRead()
            a.query_name = tag
            a.flag = record.flag
            a.query_sequence = str(Seq(duplex).reverse_complement()) if a.is_reverse else duplex
            a.reference_id = record.reference_id
            a.reference_start = record.reference_start
            a.mapping_quality = 255
            a.cigartuples = record.cigartuples
            a.next_reference_id = record.next_reference_id
            a.next_reference_start = record.next_reference_start
            a.template_length = record.template_length
            a.query_qualities = record.query_qualities
            #  Write SAM record
            out_bam_file.write(a)
    out_bam_file.close()

    # Write summary statistics
    print("Summary:")
    print("Reads processed: {}".format(reads_count))
    print("Reads discarded: {}".format(discarded_count))
    print("Reads passing filters: {}".format(reads_count_filtered))
    print("Consensus Made: {}".format(consensus_count))
    print("Duplex Made: {}".format(duplex_count))

if __name__ == "__main__":
    main()
