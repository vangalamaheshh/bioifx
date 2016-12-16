#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

__author__ = "Mahesh Vangala"
__email__  = "vangalamaheshh@gmail.com"
__date__   = "Dec, 01, 2016"

"""
Generate sequences with UMI and gene

Takes in fastq file and outputs a fastq file.
Output consists of only those sequences with
both UMI and gene sequence.

In addition, generate a map file that contains
umi and list of sequence headers
"""

import sys
import gzip
import re
from collections import defaultdict

SPACER = "GTGATTGCTTGTGACGCCTT"
GENE_DISTANCE_FROM_UMI = 21

def filter_fastq (fastq_file, spacer, gene_distance_from_umi):
    """
        Prints fastq sequences (gene sequence only)

        Trims the spacer and the umi along with oligo dTs.
        Prints fastq sequnces besides returning a dict with
        seq_to_umi info.
    """
    seq_to_umi = defaultdict(dict)
    regex_obj = re.compile(r'(?P<pre_umi>.+' + spacer + r')(?P<umi>.{10})(?P<oligo_dt>.{21})(?P<gene>.+)')
    with gzip.open(fastq_file,  'rb') as fh:
        for header in fh:
            seq = fh.readline()
            sep = fh.readline()
            ascii_line = fh.readline()
            [header, seq, sep, ascii_line] = [line.rstrip().decode('ascii') for line in
                    [header, seq, sep, ascii_line]]
            match_obj = regex_obj.search(seq)
            if match_obj:
                offset_pos = len(match_obj.group('pre_umi')) + len(match_obj.group('umi')) + \
                        len(match_obj.group('oligo_dt'))
                print(header, seq[offset_pos:], sep, ascii_line[offset_pos:], sep = "\n")
                seq_to_umi[header] = match_obj.group('umi')
    return seq_to_umi

def print_umi_info (umi_info):
    for header, umi in umi_info.items():
        sys.stderr.write(header + ',' + umi + "\n")

if __name__ == "__main__":
    umi_info = filter_fastq(sys.argv[1], SPACER, GENE_DISTANCE_FROM_UMI)
    print_umi_info(umi_info)
    sys.exit(0)
    	
