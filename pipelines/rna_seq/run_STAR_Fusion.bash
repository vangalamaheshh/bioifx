#!/bin/bash

#$ -N RUN_STAR_Fusion
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_STAR_Fusion.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_STAR_Fusion.out
#$ -pe pvm 8
#$ -l mem_free=40G 
#$ -q all.q

function map_reads {
	local input_dir=$1
	local index_dir=$2
	local not_strand_specific=$3
	local output_dir="$(dirname ${input_dir})/STAR_Fusion"
	local strand_command=""
	if [ $not_strand_specific ]; then
		echo "INFO: Processing not strand specific STAR" >&2
		strand_command="--outSAMstrandField intronMotif"
	else
		echo "INFO: Processing strand specific STAR" >&2
		strand_command="--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000"
	fi
	
	for left_mate in $(find $input_dir -mindepth 1 -maxdepth 1 -name "20*_R1*"); do
		right_mate=$(echo $left_mate | sed -e "s/_R1/_R2/")
		output_token=$(basename ${left_mate} | perl -e 'my $token = <STDIN>; chomp $token; my( $result ) = ( $token =~ /(.+?)_/ ); print $result;')
              	echo "INFO: Processing ${output_token}" >&2
		mkdir -p ${output_dir}/${output_token}
		if [ ! -f $right_mate ]; then
			local command="$STAR --runMode alignReads --runThreadN 8 --genomeDir $index_dir --readFilesIn $left_mate --readFilesCommand zcat --outFileNamePrefix ${output_dir}/${output_token}/${output_token}. --outSAMmode Full --outSAMattributes All $strand_command"	
			echo "INFO: Executing command (SE): $command" >&2
			$command
		else
			local command="$STAR --runMode alignReads --runThreadN 8 --genomeDir $index_dir --readFilesIn $left_mate $right_mate --readFilesCommand zcat --outFileNamePrefix ${output_dir}/${output_token}/${output_token}.  --outSAMmode Full --outSAMattributes All $strand_command"
			echo "INFO: Executing command (PE): $command" >&2
			$command
		fi
		sam_file="${output_dir}/${output_token}/${output_token}.Aligned.out.sam"
		bam_file=$(echo $sam_file | sed -e "s/.Aligned.out.sam/.bam/")
		sorted_bam_file=$(echo $sam_file | sed -e "s/.Aligned.out.sam/.sorted/")
		
		echo "INFO: Processing file: ${sam_file}" >&2
		$SAMTOOLS view -b -S $sam_file 1>$bam_file && $SAMTOOLS sort $bam_file $sorted_bam_file && $SAMTOOLS index ${sorted_bam_file}.bam &
		
#        	echo "INFO: Processing AddOrReplaceGroups" >&2

#       		$JAVA_EXEC -jar /zfs/cores/mbcf/mbcf-storage/devel/umv/software/picard-tools-1.115/AddOrReplaceReadGroups.jar INPUT=${sam_file} OUTPUT=${sorted_bam_file} SORT_ORDER=coordinate RGID=${output_token} RGLB=${output_token} RGPL=illumina RGPU=${output_token} RGSM=${output_token}  && echo "INFO: Processing samtools index" >&2 && $SAMTOOLS index ${sorted_bam_file} &

		# RUN STAR Fusion #
		chimeric_junction_file="${output_dir}/${output_token}/${output_token}.Chimeric.out.junction"
		chimeric_sam_file="${output_dir}/${output_token}/${output_token}.Chimeric.out.sam"
		local command="$STAR_FUSION --chimeric_junction $chimeric_junction_file -S $chimeric_sam_file --out_prefix ${output_dir}/${output_token}/${output_token}"
		echo "INFO: Executing STAR-Fusion: $command" >&2
	 	cd ${output_dir}/${output_token}/  && $command 
		for file in ${output_dir}/${output_token}/star* ; do
			out_file=$(basename $file | sed -e "s/star-fusion./${output_token}./")
			mv $file ${output_dir}/${output_token}/${out_file}
		done
	done
	
#	wait
#	mkdir -p ${output_dir}/BAM
#	for file in $(find ${output_dir} -type f -name "*sorted.ba*" ); do ln -s $file ${output_dir}/BAM/ ; done	
}

STAR="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR_2.3.0e.Linux_x86_64_static/STAR"
STAR_FUSION="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR-Fusion/STAR-Fusion"
HG19_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/STARIndex"
HG19_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
HG19_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-current/Genes/genes.gtf"
#HG19_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/STARIndex_with_ERCC"
#HG19_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_with_ERCC.fa"
#HG19_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-current/Genes/genes_with_ERCC.gtf"


MM9_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/STARIndex"
MM9_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
MM9_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-current/Genes/genes.gtf"

#PIG_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/pig/UCSC/Sus_scrofa/UCSC/susScr3/Sequence/STARIndex"
PIG_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/pig/ENSEMBL/STAR_index"
MONKEY_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/monkey/crab_eating/ENSEMBL/STAR_index"
RUN_DIR=""

JAVA_EXEC="/ifs/rcgroups/mbcf/umv/software/jdk1.8.0_45/bin/java"
SAMTOOLS="/usr/bin/samtools"


while getopts ":D:HMPKSh" opt; do
        case $opt in
                D)
			RUN_DIR=$OPTARG
			;;
		H)
			HUMAN=1
			;;
		M)
			MOUSE=1
			;;
		P)
			PIG=1
			;;
		K)	
			MONKEY=1
			;;
		S)
			NOT_STRAND_SPECIFIC=1
			;;
		h)
			echo "Usage: 	$0 -D <directory> 
				-H [this sets a flag to use human STAR index hg19]
				-P [this sets a flag to use pig STAR index]
				-M [this sets a flag to use mouse STAR index mm9]
				-K [this sets a flag to use monkey STAR index]
				-S [if you run STAR as not strand-specific; use this flag]
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

if [ -z $HUMAN ] && [ -z $MOUSE ] && [ -z $PIG ] && [ -z $MONKEY ]; then
	echo "ERROR: Either human [-H] or mouse [-M] or pig [-P] or monkey [-K] flags must be provided to the script. Exiting ..... !" >&2
	exit 1
fi

if [ $HUMAN ]; then
	echo "INFO: Loading human genome STAR index into memory" >&2
	$STAR --genomeDir $HG19_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/concat_per_sample_fastq" $HG19_INDEX $NOT_STRAND_SPECIFIC
	echo "INFO: Removing human genome STAR index from memory" >&2
	$STAR --genomeDir $HG19_INDEX --genomeLoad Remove
elif [ $MOUSE ]; then
	echo "INFO: Loading mouse genome STAR index into memory" >&2
	$STAR --genomeDir $MM9_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/concat_per_sample_fastq" $MM9_INDEX $NOT_STRAND_SPECIFIC
	echo "INFO: Removing mouse genome STAR index from memory" >&2
	$STAR --genomeDir $MM9_INDEX --genomeLoad Remove
elif [ $PIG ]; then
	echo "INFO: Loading pig genome STAR index into memory" >&2
        $STAR --genomeDir $PIG_INDEX --genomeLoad LoadAndKeep
        echo "INFO: Mapping reads" >&2
        map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/concat_per_sample_fastq" $PIG_INDEX $NOT_STRAND_SPECIFIC
        echo "INFO: Removing pig genome STAR index from memory" >&2
        $STAR --genomeDir $PIG_INDEX --genomeLoad Remove
elif [ $MONKEY ]; then
	echo "INFO: Loading monkey genome STAR index into memory" >&2
	$STAR --genomeDir $MONKEY_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/concat_per_sample_fastq" $MONKEY_INDEX $NOT_STRAND_SPECIFIC
	echo "INFO: Removing monkey genome STAR index from memory" >&2
        $STAR --genomeDir $MONKEY_INDEX --genomeLoad Remove
else
	echo "No organism selected. Exiting .... !" >&2
fi


