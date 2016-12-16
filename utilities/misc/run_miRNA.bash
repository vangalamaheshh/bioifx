#!/bin/bash

cutadapt --adapter=ATTACTCG --adapter=CGAGTAAT --minimum-length=8 --untrimmed-output=/zfs/cores/mbcf/mbcf-storage/NS500/analysis/160416_NB501431_0043_AHWKGYBGXX_Gary/analysis/cutadapt/20160417-C10ng-1-Gary_S7.notfound.fastq -o /zfs/cores/mbcf/mbcf-storage/NS500/analysis/160416_NB501431_0043_AHWKGYBGXX_Gary/analysis/cutadapt/20160417-C10ng-1-Gary_S7.clean.fastq /zfs/cores/mbcf/mbcf-storage/NS500/analysis/160416_NB501431_0043_AHWKGYBGXX_Gary/bcl2fastq/20160417-C10ng-1-Gary_S7_R1_001.fastq.gz

seqcluster collapse -f /zfs/cores/mbcf/mbcf-storage/NS500/analysis/160416_NB501431_0043_AHWKGYBGXX_Gary/analysis/cutadapt/20160417-C10ng-1-Gary_S7.clean.fastq -o /zfs/cores/mbcf/mbcf-storage/NS500/analysis/160416_NB501431_0043_AHWKGYBGXX_Gary/analysis/collapse

seqcluster seqbuster --out /zfs/cores/mbcf/mbcf-storage/NS500/analysis/160416_NB501431_0043_AHWKGYBGXX_Gary/analysis/results --hairpin /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/miRNA/hairpin.fa --mirna /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/miRNA/miRNA.str --sps hsa /zfs/cores/mbcf/mbcf-storage/NS500/analysis/160416_NB501431_0043_AHWKGYBGXX_Gary/analysis/collapse/20160417-C10ng-1-Gary_S7.clean_trimmed.fastq
