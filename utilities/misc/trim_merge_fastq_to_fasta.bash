#!/bin/bash

function filter_reads {
        local left_file=$1
	local right_file=$2
        local trimmomatic_out_dir=$3
        mkdir -p "$trimmomatic_out_dir"
        echo "Processing LEFT_MATE: ${left_file}, RIGHT_MATE: ${right_file}" >&2
                left_base=$(basename ${left_file})
                right_base=$(basename ${right_file})

                filtered_left_paired=$(echo ${trimmomatic_out_dir}/${left_base} | sed -e "s/.fastq/.trimmomatic.left_paired.fastq/")
                filtered_left_unpaired=$(echo ${trimmomatic_out_dir}/${left_base} | sed -e "s/.fastq/.trimmomatic.left_unpaired.fastq/")

                filtered_right_paired=$(echo ${trimmomatic_out_dir}/${right_base} | sed -e "s/.fastq/.trimmomatic.right_paired.fastq/")
                filtered_right_unpaired=$(echo ${trimmomatic_out_dir}/${right_base} | sed -e "s/.fastq/.trimmomatic.right_unpaired.fastq/")

                java -jar /media/mvangala/data/software/bioinfo/QC/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 \
                        ${left_file} ${right_file} \
                        $filtered_left_paired $filtered_left_unpaired \
                        $filtered_right_paired $filtered_right_unpaired \
                        ILLUMINACLIP:/media/mvangala/data/software/bioinfo/QC/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 \
                        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
		
}

function merge_reads {
        trimmomatic_dir=$1
        pandaseq_dir=$2
        mkdir -p $pandaseq_dir;
        IFS=" " read -a left_paired_files <<< $(find $trimmomatic_dir/ -type f -name "*trimmomatic.left_paired*")
        if [ $(echo "${left_paired_files[0]}" | perl -e ' my $file = <STDIN>; chomp $file; $file =~ /.gz$/ ? print STDOUT 1 : print STDOUT 0;') -eq 1 ]; then
                for left_file in ${left_paired_files[@]}; do
                        right_file=$(echo $left_file | sed -e "s/left_/right_/" | sed -e "s/_R1_/_R2_/");
                        out_token=$(basename $left_file | sed -e "s/.trimmomatic.left_paired.fastq.gz//")
                        /media/mvangala/data/software/bioinfo/format/pandaseq/pandaseq -f <(zcat $left_file) -r <(zcat $right_file) -o 6 -F -T 6 -w ${pandaseq_dir}/${out_token}.pandaseq.fastq
                done
        else
                for left_file in ${left_paired_files[@]}; do
            right_file=$(echo $left_file | sed -e "s/left_/right_/" | sed -e "s/_R1_/_R2_/");
            out_token=$(basename $left_file | sed -e "s/.trimmomatic.left_paired.fastq//")
            /media/mvangala/data/software/bioinfo/format/pandaseq/pandaseq -f $left_file -r $right_file -o 10 -F -T 6 -w ${pandaseq_dir}/${out_token}.pandaseq.fastq
        done
        fi
}

function fastq_2_fasta {
        panda_dir=$1
        fasta_dir=$2
        mkdir -p $fasta_dir
        for fastq_file in $(find $panda_dir -mindepth 1 -maxdepth 1 -type f -name "*fastq"); do
                local out_file=$(basename $fastq_file | sed -e "s/.fastq$/.fasta/")
                /media/mvangala/data/software/bioinfo/format/fastx_toolkit-0.0.14/src/fastq_to_fasta/fastq_to_fasta -i $fastq_file -o ${fasta_dir}/${out_file}
        done

}


for dir in $(find /media/mvangala/trim_and_merge/ -mindepth 1 -maxdepth 1 -type d); do
	left_file=$(find $dir -mindepth 1 -maxdepth 1 -type f -name "*_R1_*")
	right_file=$(echo $left_file | sed -e "s/_R1_/_R2_/")
	filter_reads $left_file $right_file ${dir}/trimmomatic
	merge_reads ${dir}/trimmomatic ${dir}/merge_reads_fastq
	fastq_2_fasta ${dir}/merge_reads_fastq ${dir}/merge_reads_fasta
	fastq_2_fasta ${dir} ${dir}/raw_reads_fasta
done
