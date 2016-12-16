#!/bin/bash
#$ -N FAST_QC
#$ -cwd
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/fastqc.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/fastqc.out
#$ -pe pvm 8

/zfs/cores/mbcf/mbcf-storage/devel/umv/software/FastQC/fastqc --outdir /zfs/cores/mbcf/mbcf-storage/NS500/analysis/150705_NS500305_0044_AHG5FNBGXX/fastqc/after_filtering --threads 8 /zfs/cores/mbcf/mbcf-storage/NS500/analysis/150705_NS500305_0044_AHG5FNBGXX/trimmomatic/*.fastq.gz >&2
