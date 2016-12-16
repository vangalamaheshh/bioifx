#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#-----------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Oct, 20, 2016
#-----------------------------------

import sys

def get_seq_length(fasta_file):
    """ Takes in a fasta_file (single line sequence) and returns length

        input: fasta_file <a file that contains header & sequence in each line>
        output: Returns a dict with id => seq_length

    """
    fasta_info = dict()
    with open(fasta_file, "r") as fh:
        for header in fh:
            header = header.strip()
            header = header.replace('>', '')
            seq = fh.readline().strip()
            fasta_info[header] = len(seq)
    return fasta_info

if __name__ == "__main__":
    fasta_info = get_seq_length(sys.argv[1])
    for fasta_id, seq_len in fasta_info.items():
        print(",".join([fasta_id, str(seq_len)]))
    
