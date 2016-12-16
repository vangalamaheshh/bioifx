#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#------------------------------------

samples = ["20160918_MCF7-Z_UBD3377_S26","20160918_MCF7-Q_UBD3377_S17","20160918_MCF7-R_UBD3377_S18",
           "20160918_MCF7-L_UBD3377_S12","20160918_MCF7-C_UBD3377_S3","20160918_MCF7-U_UBD3377_S21",
           "20160918_MCF7-D_UBD3377_S4","20160918_MCF7-S_UBD3377_S19","20160918_MCF7-AA_UBD3377_S27",
           #"20160918_MCF7-M_UBD3377_S13",
           "20160918_MCF7-T_UBD3377_S20","20160918_MCF7-E_UBD3377_S5",
           "20160918_MCF7-B_UBD3377_S2","20160918_MCF7-P_UBD3377_S16","20160918_MCF7-N_UBD3377_S14",
           "20160918_MCF7-F_UBD3377_S6","20160918_MCF7-V_UBD3377_S22","20160918_MCF7-A_UBD3377_S1",
           "20160918_MCF7-H_UBD3377_S8","20160918_MCF7-J_UBD3377_S10","20160918_MCF7-Y_UBD3377_S25",
           "20160918_MCF7-I_UBD3377_S9","20160918_MCF7-K_UBD3377_S11","20160918_MCF7-X_UBD3377_S24",
           "20160918_MCF7-G_UBD3377_S7","20160918_MCF7-W_UBD3377_S23","20160918_MCF7-O_UBD3377_S15"] 

rule target:
    input:
        expand( "analysis/hmm-copy/{sample}/{sample}.map.seg", sample = samples ),
        expand( "analysis/hmm-copy/{sample}/{sample}.readcounts.seg", sample = samples ),
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.gc.seg",
        expand("analysis/hmm-copy/{sample}/{sample}.corrected_CNV.pdf", sample = samples)
        

rule get_chrom_size:
    output:
        "analysis/bam2bw/hg19.Chromsizes.txt"
    params:
        "hg19"
    message: "Fetching chromosome sizes"
    shell:
        "fetchChromSizes {params} 1>{output}"


rule bam_to_bigwig:
    input:
        bam="analysis/bwa-mem/{sample}/{sample}.sorted.bam",
        chrom_size="analysis/bam2bw/hg19.Chromsizes.txt"
    output:
        protected("analysis/bam2bw/{sample}/{sample}.bw")
    params:
        "analysis/bam2bw/{sample}/{sample}"
    message: "Converting {wildcards.sample} bam to bigwig"
    shell:
        "bedtools genomecov -bg -split -ibam {input.bam} -g {input.chrom_size} 1> {params}.bg"
        " && bedSort {params}.bg {params}.sorted.bg"
        " && bedGraphToBigWig {params}.sorted.bg {input.chrom_size} {output}"

rule map_counter:
    input:
        "analysis/bam2bw/{sample}/{sample}.bw" 
    output:
        "analysis/hmm-copy/{sample}/{sample}.map.seg"
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/HMMcopy/bin/mapCounter -w 100000 "
        "{input} 1>{output}"


rule gc_counter:
    input:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
    output:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.gc.seg"
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/HMMcopy/bin/gcCounter -w 100000 "
        "{input} 1>{output}"


rule read_counter:
    input:
        "analysis/bwa-mem/{sample}/{sample}.sorted.bam"
    output:
        "analysis/hmm-copy/{sample}/{sample}.readcounts.seg"
    shell:
        "/zfs/cores/mbcf/mbcf-storage/devel/umv/software/HMMcopy/bin/readCounter -w 100000 "
        "{input} 1>{output}"

rule run_HMMcopy:
    input:
        read_file = "analysis/hmm-copy/{sample}/{sample}.readcounts.seg",
        gc_file = "/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.gc.seg",
        map_file = "analysis/hmm-copy/{sample}/{sample}.map.seg" 
    output:
        corrected_map_gc_png = "analysis/hmm-copy/{sample}/{sample}.corrected_map_gc.png",
        corrected_CNV_pdf = "analysis/hmm-copy/{sample}/{sample}.corrected_CNV.pdf"
    shell:        
        "Rscript /zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/CNV/HMMcopy.R "
        "{input.read_file} {input.gc_file} {input.map_file} {output.corrected_map_gc_png} {output.corrected_CNV_pdf}"
