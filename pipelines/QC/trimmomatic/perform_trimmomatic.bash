#!/bin/bash

#$ -N perform_trim
#$ -cwd
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/perform_trim.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/perform_trim.out
#$ -pe pvm 8

JAVA_EXEC=/ifs/rcgroups/mbcf/umv/software/jdk1.8.0_45/bin/java

function trim_reads_PE {
        local left_file=$1
        local right_file=$2
        local trim_dir=$3
	local output_token=$4
        mkdir -p $trim_dir
        filtered_left_paired="${output_token}.trimmomatic.left_paired.fastq.gz"
        filtered_left_unpaired=="${output_token}.trimmomatic.left_unpaired.fastq.gz"
        filtered_right_paired="${output_token}.trimmomatic.right_paired.fastq.gz"
        filtered_right_unpaired="${output_token}.trimmomatic.right_unpaired.fastq.gz"
        $JAVA_EXEC -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -threads 8 \
                $left_file $right_file \
                "${trim_dir}/${filtered_left_paired}" "${trim_dir}/${filtered_left_unpaired}" \
        "${trim_dir}/${filtered_right_paired}" "${trim_dir}/${filtered_right_unpaired}" \
        ILLUMINACLIP:/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
}

function trim_reads_SE {
	local left_file=$1
        local trim_dir=$2
	local output_token=$3
        mkdir -p $trim_dir
        filtered_SE_file="${output_token}.trimmomatic.single_end.fastq.gz"
        $JAVA_EXEC -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar SE -phred33 -threads 8 \
                $left_file ${trim_dir}/${filtered_SE_file} \
        ILLUMINACLIP:/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
}

output_token=$1
trim_dir=$2
left_file=$3
right_file=$4

if [ ! -f $right_file ]; then
	trim_reads_SE $left_file $trim_dir $output_token
else
	trim_reads_PE $left_file $right_file $trim_dir $output_token
fi 
