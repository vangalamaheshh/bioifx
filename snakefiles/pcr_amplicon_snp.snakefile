#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab

__author__ = "Mahesh Vangala"
__email__  = "vangalamaheshh@gmail.com"
__date__   = "Dec, 16, 2016"
__doc__    = """
    Runs bowtie and MuTect on WGS data

    
"""

import os
import sys
import pandas as pd
import snakemake
from snakemake.report import data_uri as dataURI
from snakemake.utils import report as HTMLReport

configfile: "config.yaml"
df = pd.read_csv("metasheet.csv", sep = ",", header = 0, \
        index_col = 0, comment = '#')

def getFastqFiles(wildcards):
    leftMate = df[df.index == wildcards.sample]["Filename"].tolist()[0]
    rightMate = leftMate.replace("_R1_", "_R2_")
    return ["data/" + leftMate, "data/" + rightMate ] if os.path.isfile("data/" + rightMate) \
        else ["data/" + leftMate]

def getIntervalFile(wildcards):
    intervalFile = config["params"]["mutect2"]["intervalFile"]
    return "-L " + intervalFile if intervalFile else ""

rule target:
    input:
        expand("analysis/report/alignment/{sample}/{sample}.picard.wgs_metrics.txt", \
                sample = df.index),
        "analysis/mutect2/" + config["project_name"] + ".merged.mutect2.vcf",
        "analysis/report/alignment/align_report.png"

rule runBwa:
    input:
        getFastqFiles
    output:
        protected("analysis/bwa_mem/{sample}/{sample}.sam")
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
        protected("analysis/bwa_mem/{sample}/{sample}.bam")
    message:
        "Sam to bam coversion"
    shell:
        "samtools view -bS {input} 1>{output}"


rule sortBam:
    input:
        unsortedBam = "analysis/bwa_mem/{sample}/{sample}.bam"
    output:
        sortedBam = protected("analysis/bwa_mem/{sample}/{sample}.sorted.bam"),
        bam_index = protected("analysis/bwa_mem/{sample}/{sample}.sorted.bam.bai")
    message:
        "Sorting and indexing bam"
    shell:
        "samtools sort -f {input.unsortedBam} {output.sortedBam} "
        "&& samtools index {output.sortedBam}"

rule samtoolsStats:
    input:
        unsortedBam = "analysis/bwa_mem/{sample}/{sample}.bam"
    output:
        samStats = protected("analysis/report/alignment/{sample}/{sample}.samtools.stats.txt")
    message:
        "Running samtools stats on {wildcards.sample}"
    shell:
        "samtools stats {input.unsortedBam} | grep ^SN | "
        "gawk 'BEGIN {{ FS=\"\t\"; }} {{ print $2,$3; }}' 1>{output.samStats}"

rule picardStats:
    input:
        sortedBam = "analysis/bwa_mem/{sample}/{sample}.sorted.bam"
    output:
        picardStats = protected("analysis/report/alignment/{sample}/" + \
                        "{sample}.picard.wgs_metrics.txt")
    message:
        "Running picard wgs stats on {wildcards.sample}"
    params:
        refFasta = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/ref_files/" \
                    + "human/ucsc/hg19/Sequence/WholeGenomeFasta/genome.fa"
    shell:
        "picard CollectWgsMetrics I={input.sortedBam} O={output.picardStats} "
        "R={params.refFasta}" 

rule mapReportMatrix:
    input:
        metricsList = expand("analysis/report/alignment/{sample}/" + \
                "{sample}.samtools.stats.txt", sample = df.index)
    output:
        csv = "analysis/report/alignment/align_report.csv"
    message:
        "Gather samtools stats into csv"
    run:
        argList = " -s " + " -s ".join(input.metricsList)
        shell("perl /ifs/rcgroups/ccgd/ccgd-data/home/umv/git/bioifx/utilities/alignment"
        + "/sam_stats_matrix.pl " + argList + " 1>{output.csv}") 
        
rule mapReportPlot:
    input:
        csv = "analysis/report/alignment/align_report.csv"
    output:
        png = "analysis/report/alignment/align_report.png"
    message:
        "Plotting alignment PNG"
    shell:
        "Rscript /ifs/rcgroups/ccgd/ccgd-data/home/umv/git/bioifx/utilities/"
        "alignment/sam_stats_matrix.R {input.csv} {output.png}"

rule runMuTect2:
    input:
        "analysis/bwa_mem/{sample}/{sample}.sorted.bam"
    output:
        protected("analysis/mutect2/{sample}/{sample}.mutect2.vcf")
    message:
        "Running MuTect2 on {wildcards.sample}"
    threads: 8
    params:
        gatkexec = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/software/GenomeAnalysisTK.jar",
        refFasta = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/ref_files/" \
                    + "human/ucsc/hg19/Sequence/WholeGenomeFasta/genome.fa",
        dbSnp = "/ifs/rcgroups/ccgd/ccgd-data/home/umv/software/" \
                    + "muTect-1.1.4/lib/dbsnp_138.hg19.vcf",
        intervalFile = getIntervalFile
    shell:
        "java -jar {params.gatkexec} -T MuTect2 -nct {threads} -R {params.refFasta} "
        "-I:tumor {input} --dbsnp {params.dbSnp} {params.intervalFile} -o {output}"

rule mergeVCFs:
    input:
        vcfList = expand("analysis/mutect2/{sample}/{sample}.mutect2.vcf", sample = df.index)
    output:
        mergedVcf = "analysis/mutect2/" + config["project_name"] + ".merged.mutect2.vcf"
    message:
        "Merging vcfs into one"
    threads: 4
    run:
        vcfs = " --variant " + " --variant ".join(input.vcfList)
        shell("java -jar /ifs/rcgroups/ccgd/ccgd-data/home/umv/software/GenomeAnalysisTK.jar " +
            "-T CombineVariants -R /ifs/rcgroups/ccgd/ccgd-data/home/umv/ref_files/human/ucsc/" +
            "hg19/Sequence/WholeGenomeFasta/genome.fa " + vcfs + " -o " +
            output.mergedVcf + " -genotypeMergeOptions UNIQUIFY -nt 4")

rule TiTvRatio:
    input:
        csv = "analysis/report/TiTv/{sample}/csv/{sample}.ti_tv_ratio.csv"
    output:
        png = "analysis/report/TiTv/{sample}/png/{sample}.ti_tv_ratio.png"
    message:
        "Plotting Transitions to Transversions Ratio"
    shell:
        "Rscript /ifs/rcgroups/ccgd/ccgd-data/home/umv/git/bioifx/utilities/SNP/ti_tv_ratio.R "
        "{input.csv} {output.png}"

rule generateReport:
    output:
        html = "analysis/report/" + config["project_name"] + ".html"
    message:
        "Generating HTML Report"
    run:
        report = """
=====================
PCR Amplicon Project
=====================

Workflow:
=========
Alignment:
==========
    Raw reads are mapped to human reference genome (UCSC build - hg19) using *BWA mem* aligner. The alignment statistics are generated using *samtools stats* and *picard wgs_metrics*. The alignment stats are given in the report below. 
""" 

        report += "\n\t.. image:: " + dataURI("analysis/report/" +
                "alignment/align_report.png") + "\n" 


        report += """
SNP calling:
============
    Aligned bam files are processing using MuTect2 software (now part of GATK). The resultant vcf files are merged into one vcf file for archival purposes. For further analysis, SNPs are filtered using *PASS* criterion for each sample. Transitions to Transversions ratio is generated for every 100,000bp window for *chr7, chr17 and chrX*. 
"""
        report += "\n\t.. image:: " + dataURI("analysis/report/ti_tv_ratio.png") + "\n"

        HTMLReport(report, output.html, 
            metadata = "Center for Cancer Genome Discovery",
            **{'Copyrights': "/ifs/rcgroups/ccgd/ccgd-data/home/umv/misc/ccgd.jpg"})
