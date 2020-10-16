#!/usr/bin/env python
"""
This tool creates consensus reads using aligned paired-end reads generated with
the Duplex Sequencing method. The UMIs must be in the headers like this:

HWI-ST1390:273:CDTB8ACXX:6:2316:21033:100486_ATACTGTCAATT
HWI-ST1390:273:CDTB8ACXX:6:2316:21033:100486_ATACTGTCAATT
...

This is a simplified, modified version of the code present in:

https://github.com/Kennedy-Lab-UW/Duplex-Sequencing

Input:
	A position-sorted paired-end BAM file containing reads with a duplex UMI in the header.

Output:
	1: A pair of FASTQ files containing the paired duplexs
	2: A pair of FASTQ files containing the unpaired duplexs
	3: A pair of FASTQ files containing the consensus reads
	4: A BAM file containing discarded reads

Author: jc.fernandez.navarro@gmail.com
"""
import pysam
import gzip
from fastqutils import *
from collections import defaultdict, Counter
from argparse import ArgumentParser, RawDescriptionHelpFormatter


# TODO add option to clusters tags allowing for mismatches
# TODO add option to cluster tags with a window size in position

def consensus_maker(grouped_reads_list, cut_off):
    """
    Create a consensus read using the most_common letter from each position
    """
    num_reads = float(len(grouped_reads_list))
    counters = [Counter(x) for x in zip(*grouped_reads_list)]
    consensus_read = ''.join(
        [c.most_common()[0][0] if (c.most_common()[0][1] / num_reads) >= cut_off else 'N' for c in counters])
    return consensus_read


def consensus_quality(qual_list):
    """
    Create a consensus quality using the average of the qualities
    """
    qual_avg = [int(sum(qual_score) / len(qual_list)) for qual_score in zip(*qual_list)]
    return [x if x < 41 else 41 for x in qual_avg]


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('infile', type=str, help="input BAM file")
    parser.add_argument('--min-reads', type=int, default=2, dest='min_reads',
                        help="Minimum number of reads allowed to comprise a consensus. [2]")
    parser.add_argument('--max-reads', type=int, default=1000, dest='max_reads',
                        help="Maximum number of reads allowed to comprise a consensus. [1000]")
    parser.add_argument('--homo-count', type=int, default=8, dest='rep_filt',
                        help="Maximum number of homopolymer supported in the UMI [A, C, G, T]. [8]")
    parser.add_argument('--cutoff', type=float, default=.7, dest='cut_off',
                        help="Percentage of nucleotides at a given position in a read that must be identical in order "
                             "for a consensus to be called at that position. [0.7]")
    parser.add_argument('--filter-soft-clip', action="store_true", default=False, dest='softclip_filter',
                        help="Discard reads that are soft-clipped in their alignment.")
    parser.add_argument('--filter-overlap', action="store_true", default=False, dest='overlap_filter',
                        help="Discard reads overlapping with their pair-end mate.")
    parser.add_argument('--filter-pair', action="store_true", default=False, dest='pair_filter',
                        help="Discard reads that do not have a proper pair-end mate.")
    parser.add_argument('--filter-secondary', action="store_true", default=False, dest='secondary_filter',
                        help="Discard reads that are not primary alignment.")
    parser.add_argument('--filter-singleton', action="store_true", default=False, dest='single_filter',
                        help="Discard reads that are singletons (only one pair is aligned).")
    parser.add_argument('--disable-pos-clustering', action="store_true", default=False, dest='no_pos_clustering',
                        help="Do not take the read start position when clustering UMIs.")
    parser.add_argument('--maxN', type=float, default=.5, dest='Ncut_off',
                        help="Maximum fraction of Ns allowed in a consensus [0.5]")
    o = parser.parse_args()

    # Variables
    reads_count = 0
    reads_count_filtered = 0
    consensus_count = 0
    duplex_count = 0
    simplex_count = 0
    read_dict = defaultdict(lambda: defaultdict(list))
    min_reads = o.min_reads
    max_reads = o.max_reads
    homo_count = o.rep_filt
    cut_off = o.cut_off
    Ncut_off = o.Ncut_off
    do_overlap_filter = o.overlap_filter
    do_softclip_filter = o.softclip_filter
    do_secondary = o.secondary_filter
    do_pair_filter = o.pair_filter
    do_single_filter = o.single_filter
    do_N_filter = Ncut_off <= 1.0 and Ncut_off >= 0.0

    # Start going through the input BAM file, one position at a time to group reads by start_pos + tag
    in_bam_file = pysam.Samfile(o.infile, "rb")
    discarded_bam_file = pysam.Samfile('discarded.bam', "wb", template=in_bam_file)
    print("Grouping reads by position-tag...")
    for record in in_bam_file.fetch(until_eof=True):
        reads_count += 1

        #  Get tag
        tag = record.query_name.split('_')[-1]
        tag += ":1" if record.is_read1 else ":2"

        # Mapping filters
        unmapped = record.is_unmapped
        singleton = do_single_filter and record.mate_is_unmapped
        not_proper_pair = do_pair_filter and not record.is_proper_pair
        secondary = do_secondary and record.is_secondary

        # Overlap filter
        overlap = False
        if do_overlap_filter and not not_proper_pair and not singleton:
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
        if any([overlap, soft_clip, homo_A, homo_C, homo_G, homo_T, singleton, not_proper_pair, unmapped, secondary]):
            discarded_bam_file.write(record)
        else:
            reads_count_filtered += 1
            key_pos = 0 if o.no_pos_clustering else record.reference_start
            read_dict[key_pos][tag].append(record)

    in_bam_file.close()
    discarded_bam_file.close()

    # Iterate grouped reads to build consensus reads
    print("Creating consensus reads...")
    out_R1_handle = gzip.open("duplex_R1.fq.gz", 'wt')
    out_R1_writer = writefq(out_R1_handle)
    out_R1_consensus_handle = gzip.open("consensus_R1.fq.gz", 'wt')
    out_R1_consensus_writer = writefq(out_R1_consensus_handle)
    out_R1_simplex_handle = gzip.open("simplex_R1.fq.gz", 'wt')
    out_R1_simplex_writer = writefq(out_R1_simplex_handle)
    out_R2_handle = gzip.open("duplex_R2.fq.gz", 'wt')
    out_R2_writer = writefq(out_R2_handle)
    out_R2_consensus_handle = gzip.open("consensus_R2.fq.gz", 'wt')
    out_R2_consensus_writer = writefq(out_R2_consensus_handle)
    out_R2_simplex_handle = gzip.open("simplex_R2.fq.gz", 'wt')
    out_R2_simplex_writer = writefq(out_R2_simplex_handle)
    for pos, tag_records in read_dict.items():
        consensus_dict = defaultdict(lambda: defaultdict(list))
        for tag, records in tag_records.items():
            sequences = [x.query_sequence for x in records]
            qualities = [x.query_qualities for x in records]
            num_sequences = len(sequences)
            first_length = len(sequences[0])
            if num_sequences >= min_reads and num_sequences <= max_reads and \
                    all(len(x) == first_length for x in sequences):
                consensus = consensus_maker(sequences, cut_off)
                consensus_qual = consensus_quality(qualities)
                # Filter out consensuses with too many Ns in them
                if do_N_filter and (consensus.count("N") / float(len(consensus)) >= Ncut_off):
                    continue
                consensus_count += 1
                #  Store records (there must be one record per tag)
                clean_tag, info = tag.split(":")
                consensus_dict[clean_tag][info] = (consensus, consensus_qual)
                if info == '1':
                    out_R1_consensus_writer.send(("{}:{}/1".format(clean_tag, num_sequences),
                                                  consensus, ''.join(chr(x + 33) for x in consensus_qual)))
                else:
                    out_R2_consensus_writer.send(("{}:{}/2".format(clean_tag, num_sequences),
                                                  consensus, ''.join(chr(x + 33) for x in consensus_qual)))
        # Iterate consensus dict to find duplexs otherwise write simplex
        processed = set()
        for tag, record in consensus_dict.items():
            #  Create reversed tag
            switch_tag = "{}{}".format(tag[int(len(tag) / 2):],
                                       tag[:int(len(tag) / 2)])

            duplex_read1 = None
            duplex_read2 = None
            if switch_tag in consensus_dict and switch_tag not in processed:
                if len(record['1']) != 0 and len(consensus_dict[switch_tag]['2']) != 0:
                    duplex_seq = consensus_maker([record['1'][0],
                                                  consensus_dict[switch_tag]['2'][0]],
                                                 1.0)
                    duplex_qual = consensus_quality([record['1'][1],
                                                     consensus_dict[switch_tag]['2'][1]])
                    # Filter out duplex with too many Ns in them
                    if not (do_N_filter and (duplex_seq.count("N") / float(len(duplex_seq)) >= Ncut_off)):
                        duplex_read1 = (tag + "/1", duplex_seq, ''.join(chr(x + 33) for x in duplex_qual))

                if len(record['2']) != 0 and len(consensus_dict[switch_tag]['1']) != 0:
                    duplex_seq = consensus_maker([record['2'][0],
                                                  consensus_dict[switch_tag]['1'][0]],
                                                 1.0)
                    duplex_qual = consensus_quality([record['2'][1],
                                                     consensus_dict[switch_tag]['1'][1]])
                    # Filter out duplex with too many Ns in them
                    if not (do_N_filter and (duplex_seq.count("N") / float(len(duplex_seq)) >= Ncut_off)):
                        duplex_read2 = (tag + "/2", duplex_seq, ''.join(chr(x + 33) for x in duplex_qual))

                #  To not visit this duplex again
                processed.add(switch_tag)
                processed.add(tag)

            if duplex_read1 is not None and duplex_read2 is not None:
                duplex_count += 1
                out_R1_writer.send(duplex_read1)
                out_R2_writer.send(duplex_read2)
            elif duplex_read1 is not None:
                simplex_count += 1
                out_R1_simplex_writer.send(duplex_read1)
            elif duplex_read2 is not None:
                simplex_count += 1
                out_R2_simplex_writer.send(duplex_read2)

    out_R1_handle.close()
    out_R2_handle.close()
    out_R1_consensus_handle.close()
    out_R2_consensus_handle.close()
    out_R1_simplex_handle.close()
    out_R2_simplex_handle.close()

    # Write summary statistics
    print("Summary:")
    print("Reads processed: {}".format(reads_count))
    print("Reads discarded: {}".format(reads_count - reads_count_filtered))
    print("Reads passing filters: {}".format(reads_count_filtered))
    print("Consensus made (paired and unpaired): {}".format(consensus_count))
    print("Duplex paired made: {}".format(duplex_count))
    print("Duplex unpaired made: {}".format(simplex_count))


if __name__ == "__main__":
    main()
