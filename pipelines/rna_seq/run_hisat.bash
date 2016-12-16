#!/bin/bash

#$ -N RUN_HISAT
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_hisat.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_hisat.out
#$ -pe pvm 10
#$ -l mem_free=30G
#$ -q all.q

function map_reads {
        local input_dir=$1
        local index_info=$2
	local not_strand_specific=$3
	local output_dir="${input_dir}/HISAT_with_ERCC"
	
	if [ $not_strand_specific ]; then
                echo "INFO: Processing not strand specific HISAT" >&2
                strand_command="--rna-strandness unstranded"
        else
                echo "INFO: Processing strand specific HISAT" >&2
                strand_command="--rna-strandness R"
        fi


	for left_mate in $(find ${input_dir}/concat_per_sample_fastq -mindepth 1 -maxdepth 1 -name "[0-9]*_R1*fastq.gz"); do
                right_mate=$(echo $left_mate | sed -e "s/_R1/_R2/")
		output_token=$(basename ${left_mate} | sed -e "s/_.*//")
		echo "INFO: Processing ${output_token}" >&2
                mkdir -p ${output_dir}/${output_token}
		out_prefix="${output_dir}/${output_token}/${output_token}"
                if [ ! -f $right_mate ]; then
                        local command="$HISAT_EXEC --rg-id ${output_token} --rg PL:illumina --rg LB:${output_token} --rg SM:${output_token} $strand_command --threads 10 -x ${index_info} -U $left_mate -S ${out_prefix}.sam 2>${out_prefix}.log"
                        echo "INFO: Executing command (SE): $command" >&2
                        $command
                else
                        local command="$HISAT_EXEC --rg-id ${output_token} --rg PL:illumina --rg LB:${output_token} --rg SM:${output_token} $strand_command --threads 10 -x ${index_info} -1 $left_mate -2 $right_mate -S ${out_prefix}.sam 2>${out_prefix}.log"
                        echo "INFO: Executing command (PE): $command" >&2
                        $command
                fi
                sam_file="${out_prefix}.sam"
                bam_file=$(echo $sam_file | sed -e "s/.sam/.bam/")
		sorted_bam_file=$(echo $sam_file | sed -e "s/.sam/.sorted/")
		/usr/bin/samtools view -b -S $sam_file 1>$bam_file && /usr/bin/samtools sort $bam_file $sorted_bam_file && /usr/bin/samtools index ${sorted_bam_file}.bam &
	done
	wait
#        mkdir -p ${output_dir}/BAM
#        for file in $(find ${output_dir} -type f -name "*sorted.ba*" ); do ln -s $file ${output_dir}/BAM/ ; done
}

HISAT_EXEC="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/hisat-0.1.6-beta/hisat"

MM9_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/HISATIndex/mm9_hisat"
MM9_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
MM9_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-current/Genes/genes.gtf"

#HG19_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/HISATIndex/hg19_hisat"
HG19_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/HISATIndex_with_ERCC/hg19_with_ERCC"

while getopts ":D:HMSh" opt; do
        case $opt in
                D)
                        RUN_DIR=$OPTARG
                        ;;
                M)
                        MOUSE=1
                        ;;
		H)
			HUMAN=1
			;;
		S)
                        NOT_STRAND_SPECIFIC=1
                        ;;
                h)
                        echo "Usage:    $0 
				-D </absolute/path/to/run_dir>
				-H [this sets a flag to use human HISAT index hg19]
                                -M [this sets a flag to use mouse HISAT index mm9]
				-S [if you want to run HISAT as not strand-specific; use this flag]
                                -h [prints this help message and exits]"
                        exit 1
                        ;;
                \?)
                        echo "Invalid option: $OPTARG"
                        exit 1
                        ;;
                :)
                        echo "$OPTARG - Requires argument"
                        exit 1
                        ;;
        esac
done

if [ $MOUSE ]; then
	echo "INFO: Processing mouse mapping" >&2
	map_reads $RUN_DIR $MM9_INDEX $NOT_STRAND_SPECIFIC
elif [ $HUMAN ]; then
	echo "INFO: Processing human mapping" >&2
	map_reads $RUN_DIR $HG19_INDEX $NOT_STRAND_SPECIFIC
fi
