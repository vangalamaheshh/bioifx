#!/usr/bin/env python
#vim: syntax=python tabstop=2 expandtab

import snakemake
import pandas as pd
import sys

meta = pd.read_csv("metasheet.csv", header = 0, sep = ",", index_col = 0)

def getFullPath(wildcards):
  return meta[meta.index.values == wildcards.sample]["FileName"]

rule all:
  input:
    expand("analysis/fastq/{sample}/{sample}.left.fastq", sample = meta.index.values)

rule bam2fastq:
  input:
    bamFile = getFullPath
  output:
    fastqLeftFile = "analysis/fastq/{sample}/{sample}.left.fastq",
    fastqRightFile = "analysis/fastq/{sample}/{sample}.right.fastq"
  message:
    "Converting PE Bam file to Fastq"
  shell:
    "picard SamToFastq I={input.bamFile} F={output.fastqLeftFile} "
    "F2={output.fastqRightFile} "
