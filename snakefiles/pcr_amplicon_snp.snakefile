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
        bwaexec = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/software/bwa-0.7.12/bwa",
        sampleName = lambda wildcards: wildcards.sample
    shell:
        "{params.bwaexec} mem -M -t {threads} "
        "-R '@RG\\tID:{params.sampleName}\\tLB:{params.sampleName}\\tPU:Illumina\\tSM:{params.sampleName}' "
        "{params.bwaIndex} {input} 1>{output}"
       
rule samToBam:
    input:
        "analysis/bwa_mem/{sample}/{sample}.sam"
    output:
        "analysis/bwa_mem/{sample}/{sample}.bam"
    message:
        "Sam to bam coversion"
    shell:
        "samtools view -bS {input} 1>{output}"


rule sortBam:
    input:
        unsortedBam = "analysis/bwa_mem/{sample}/{sample}.bam"
    output:
        sortedBam = "analysis/bwa_mem/{sample}/{sample}.sorted.bam",
        bam_index = "analysis/bwa_mem/{sample}/{sample}.sorted.bam.bai"
    message:
        "Sorting and indexing bam"
    shell:
        "samtools sort -f {input.unsortedBam} {output.sortedBam} "
        "&& samtools index {output.sortedBam}"

rule runMuTect2:
    input:
        "analysis/bwa_mem/{sample}/{sample}.sorted.bam"
    output:
        "analysis/mutect2/{sample}/{sample}.mutect2.vcf"
    message:
        "Running MuTect2 on {wildcards.sample}"
    threads: 8
    params:
        gatkexec = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/software/GenomeAnalysisTK.jar",
        refFasta = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/ref_files/" \
                    + "human/ucsc/hg19/Sequence/WholeGenomeFasta/genome.fa",
        dbSnp = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/software/" \
                    + "muTect-1.1.4/lib/dbsnp_138.hg19.vcf"
    shell:
        "java -jar {params.gatkexec} -T MuTect2 -nct {threads} -R {params.refFasta} "
        "-I:tumor {input} --dbsnp {params.dbSnp} -o {output}"

 
