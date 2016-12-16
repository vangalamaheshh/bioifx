#!/bin/bash

#$ -cwd
#$ -N BWA_MEM
#$ -e bwa_mem.log
#$ -o bwa_mem.out
#$ -pe pvm 8

NS500_DIR=/zfs/cores/mbcf/mbcf-storage/NS500/analysis/160918_NB501431_0164_AHLG3CBGXY_UBD3377
OUT_DIR="${NS500_DIR}/analysis/bwa-mem"

mkdir -p $OUT_DIR

for left_mate in $(find ${NS500_DIR}/data/ -type l -name "[0-9]*_R1_*fastq.gz"); do
	echo "INFO: Processing $left_mate" >&2
	echo "-----------------------" >&2
	right_mate=$(basename $left_mate | sed -e "s/_R1_/_R2_/")
	right_mate=${NS500_DIR}/data/${right_mate}
	out_token=$(basename $left_mate | sed -e "s/_R.*_001.fastq.gz//")
	mkdir -p ${OUT_DIR}/${out_token}
	#/zfs/cores/mbcf/mbcf-storage/devel/umv/software/bwa-0.7.12/bwa mem -M -t 8 /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa $left_mate $right_mate 1>${OUT_DIR}/${out_token}/${out_token}.sam
	#samtools view -b -S ${OUT_DIR}/${out_token}/${out_token}.sam  1>${OUT_DIR}/${out_token}/${out_token}.bam && \
	#samtools sort  ${OUT_DIR}/${out_token}/${out_token}.bam 1> ${OUT_DIR}/${out_token}/${out_token}.sorted.bam && \
	#samtools index ${OUT_DIR}/${out_token}/${out_token}.sorted.bam && \
	#samtools stats ${OUT_DIR}/${out_token}/${out_token}.bam | grep ^SN | gawk 'BEGIN{ FS="\t"; }{ print $2,$3; }' 1>${OUT_DIR}/${out_token}/${out_token}.stats.txt &
	#------------------------------
	# Picard stats on whole genome
	#------------------------------
	picard CollectWgsMetrics I=${OUT_DIR}/${out_token}/${out_token}.sorted.bam O=${OUT_DIR}/${out_token}/${out_token}.picard.wgs_metrics.txt R=/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa &
done

wait
echo "Done!" >&2
