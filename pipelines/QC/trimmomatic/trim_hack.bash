#!/bin/bash

#$ -N perform_trim
#$ -cwd
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/perform_trim.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/perform_trim.out
#$ -pe pvm 8

ANALYSIS_DIR=/zfs/cores/mbcf/mbcf-storage/NS500/analysis/150605_NS500305_0022_AHGC27BGXX

function trim_reads {
        local left_file=$1
        local right_file=$2
        local trim_dir=$3
        mkdir -p $trim_dir
	trim_log=$(basename $left_file | sed -e "s/_.*//");
        filtered_left_paired=$(basename $left_file | sed -e "s/_R1.fastq/.trimmomatic.left_paired.fastq/")
        filtered_left_unpaired==$(basename $left_file | sed -e "s/_R1.fastq/.trimmomatic.left_unpaired.fastq/")
        filtered_right_paired=$(basename $right_file | sed -e "s/_R2.fastq/.trimmomatic.right_paired.fastq/")
        filtered_right_unpaired=$(basename $right_file | sed -e "s/_R2.fastq/.trimmomatic.right_unpaired.fastq/")
        java -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -threads 8 -trimlog "${trim_dir}/${trim_log}.trim.log"\
                $left_file $right_file \
                "${trim_dir}/${filtered_left_paired}" "${trim_dir}/${filtered_left_unpaired}" \
        "${trim_dir}/${filtered_right_paired}" "${trim_dir}/${filtered_right_unpaired}" \
        ILLUMINACLIP:/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
}


for left_file in $(find /zfs/cores/mbcf/mbcf-storage/NS500/analysis/150605_NS500305_0022_AHGC27BGXX/concat_hack/ -mindepth 1 -maxdepth 1 -type f -name "20*R1*fastq.gz"); do
        echo "INFO:Processing sample: $left_file" >&2
        echo "-----------------------------------" >&2
        right_file=$(echo $left_file | sed -e "s/_R1/_R2/")
        sample_name=$(basename $left_file | sed -e "s/_R1.fastq.gz//");
        trim_reads $left_file $right_file "${ANALYSIS_DIR}/trimmomatic_hack/${sample_name}"
done
