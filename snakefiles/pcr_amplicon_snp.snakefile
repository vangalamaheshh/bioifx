#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

__author__ = "Mahesh Vangala"
__email__  = "vangalamaheshh@gmail.com"
__date__   = "Dec, 16, 2016"
__doc__    = """
    Runs bowtie and MuTect on WGS data

    
"""

import os
import pandas as pd
import snakemake

df = pd.read_csv("metasheet.csv", sep = ",", header = 0, \
        index_col = 0, comment = '#')

def getFastqFiles(wildcards):
    leftMate = df[df.index == wildcards.sample]["Filename"].tolist()[0]
    rightMate = leftMate.replace("_R1_", "_R2_")
    return ["data/" + leftMate, "data/" + rightMate ] if os.path.isfile("data/" + rightMate) \
        else ["data/" + leftMate]

rule target:
    input:
        expand("analysis/bwa_mem/{sample}/{sample}.sam", sample = df.index)

rule runBwa:
    input:
        getFastqFiles
    output:
        "analysis/bwa_mem/{sample}/{sample}.sam"
    threads: 8
    message:
        "Running Bwa on {wildcards.sample}"
    params:
        bwaIndex = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/" + \
                        "ref_files/human/ucsc/hg19/Sequence/BWAIndexBuilt/hg19",
        bwaexec = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/software/bwa-0.7.12/bwa"
    shell:
        "{params.bwaexec} mem -M -t {threads} {params.bwaIndex} {input} 1>{output}"
        
