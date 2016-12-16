#!/bin/bash

#$ -N PandaSeq
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/pandaseq/pandaseq.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/pandaseq/pandaseq.out

in_dir="/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/140924_M00851_0283_000000000-AALLB/demultiplex_run"
out_dir="/zfs/cores/mbcf/mbcf-storage/devel/umv/pandaseq/kai"

for left_mate_file in $(find $in_dir -type f -name "*_R1*.fastq.gz"); do
	right_mate_file=$(echo $left_mate_file | sed -e "s/_R1_/_R2_/")
	token=$(basename $left_mate_file | sed -e "s/_L001_R1_001.fastq.gz//")
	echo "INFO: Processing $left_mate_file" >&2
	
	/usr/local/bin/pandaseq -f $left_mate_file -r $right_mate_file -6 -F -l 25 -T 8 1>${out_dir}/${token}.pandaseq.fastq
done
