#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

#----------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Nov, 1st, 2016
#-----------------------------------------

import sys, re
import subprocess
import io

def get_trim_position(seq, token="AAAAAAAAAA"):
    trim_pos = -1
    regex_obj = re.compile(token)
    match_obj = regex_obj.search(seq)
    if match_obj:
        trim_pos = match_obj.start() 
    return trim_pos


def print_trimmed_reads():
    for line in sys.stdin:
        header = line.strip()
        seq = sys.stdin.readline().strip()
        sep = sys.stdin.readline().strip()
        ascii_info = sys.stdin.readline().strip()
        trim_pos = get_trim_position(seq)
        if trim_pos < 0:
            trim_pos = len(seq)
        if trim_pos >= 18:
            print(header, seq[:trim_pos], sep, ascii_info[:trim_pos], sep = "\n") 
                

if __name__ == "__main__":
    print_trimmed_reads()

