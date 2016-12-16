#!/bin/bash

#$ -N sort_bam
#$ -cwd
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/sort_bam.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/sort_bam.out


#for token in mm9; do
	for bam_file in $(find ${1} -type f -name "*.bam"); do 
		echo "Processing $bam_file";
		#bam_file=$(echo $file | sed -e "s/sam/bam/")
		sorted_file=$(echo $bam_file | sed -e "s/bam/sorted/")
		#samtools view -b -S $file 1>$bam_file
		samtools sort  $bam_file 1>${sorted_file}.bam
		samtools index ${sorted_file}.bam
	done
#done
