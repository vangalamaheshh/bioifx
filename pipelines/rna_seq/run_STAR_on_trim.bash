#!/bin/bash

#$ -N RUN_STAR
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_STAR.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_STAR.out
#$ -pe pvm 8
#$ -l mem_free=40G 
#$ -q all.q

function map_reads {
	local input_dir=$1
	local index_dir=$2
	local not_strand_specific=$3
	local output_dir="$(dirname ${input_dir})/STAR_with_TRIM"
	local strand_command=""
	if [ $not_strand_specific ]; then
		echo "INFO: Processing not strand specific STAR" >&2
		strand_command="--outSAMstrandField intronMotif"
	else
		echo "INFO: Processing strand specific STAR" >&2
		strand_command="--outFilterIntronMotifs RemoveNoncanonical"
	fi
	
	for left_mate in $(find $input_dir -mindepth 1 -maxdepth 1 -name "*single_end.fastq.gz" -o -name "*left_paired.fastq.gz"); do
		right_mate=$(echo $left_mate | sed -e "s/left_paired/right_paired/")
		output_token=$(basename ${left_mate} | perl -e 'my $token = <STDIN>; chomp $token; my( $result ) = ( $token =~ /(.+?)_/ ); print $result;')
              	echo "INFO: Processing ${output_token}" >&2
		mkdir -p ${output_dir}/${output_token}
		if [ ! -f $right_mate ]; then
			local command="$STAR --runMode alignReads --runThreadN 8 --genomeDir $index_dir --readFilesIn $left_mate --readFilesCommand zcat --outFileNamePrefix ${output_dir}/${output_token}/${output_token}. --outSAMmode Full --outSAMattributes All  ${strand_command} --outSAMattrRGline ID:${output_token} PL:illumina LB:${output_token} SM:${output_token} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"	
			echo "INFO: Executing command (SE): $command" >&2
			$command
		else
			local command="$STAR --runMode alignReads --runThreadN 8 --genomeDir $index_dir --readFilesIn $left_mate $right_mate --readFilesCommand zcat --outFileNamePrefix ${output_dir}/${output_token}/${output_token}. --outSAMmode Full --outSAMattributes All $strand_command --outSAMattrRGline ID:${output_token} PL:illumina LB:${output_token} SM:${output_token} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"
			echo "INFO: Executing command (PE): $command" >&2
			$command
		fi
		mv "${output_dir}/${output_token}/${output_token}.Aligned.sortedByCoord.out.bam" "${output_dir}/${output_token}/${output_token}.sorted.bam"
		mv "${output_dir}/${output_token}/${output_token}.ReadsPerGene.out.tab" "${output_dir}/${output_token}/${output_token}.counts.tab"
		$SAMTOOLS index "${output_dir}/${output_token}/${output_token}.sorted.bam" &	
#		$SAMTOOLS view -b -S $sam_file 1>$bam_file && $SAMTOOLS sort $bam_file $sorted_bam_file && $SAMTOOLS index ${sorted_bam_file}.bam &
	done
	
	wait
	mkdir -p ${output_dir}/BAM
	cd ${output_dir}/BAM
	for file in $(find ../$(basename ${output_dir})/ -type f -name "*sorted.ba*" ); do ln -s $file ${output_dir}/BAM/ ; done	
	STAR_LOG_FILES=""
	for file in $(find ${output_dir} -type f -name "*.Log.final.out"); do STAR_LOG_FILES="${STAR_LOG_FILES} -l $file "; done
	perl $STAR_REPORT $STAR_LOG_FILES 1>${output_dir}/STAR_Report.csv
}

STAR="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR"
STAR_REPORT="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/report_generation/STAR_reports.pl"

#HG19_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/STARIndex"
#HG19_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
#HG19_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-current/Genes/genes.gtf"
HG19_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/STARIndex_with_ERCC"
HG19_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_with_ERCC.fa"
HG19_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-current/Genes/genes_with_ERCC.gtf"


#MM9_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/STARIndex"
#MM9_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
#MM9_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-current/Genes/genes.gtf"
MM9_INDEX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/STARIndex_with_ERCC"
MM9_GENOME="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome_with_ERCC.fa"
MM9_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-current/Genes/genes_with_ERCC.gtf"


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
	map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/trimmomatic" $HG19_INDEX $NOT_STRAND_SPECIFIC
	echo "INFO: Removing human genome STAR index from memory" >&2
	$STAR --genomeDir $HG19_INDEX --genomeLoad Remove
elif [ $MOUSE ]; then
	echo "INFO: Loading mouse genome STAR index into memory" >&2
	$STAR --genomeDir $MM9_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/trimmomatic" $MM9_INDEX $NOT_STRAND_SPECIFIC
	echo "INFO: Removing mouse genome STAR index from memory" >&2
	$STAR --genomeDir $MM9_INDEX --genomeLoad Remove
elif [ $PIG ]; then
	echo "INFO: Loading pig genome STAR index into memory" >&2
        $STAR --genomeDir $PIG_INDEX --genomeLoad LoadAndKeep
        echo "INFO: Mapping reads" >&2
        map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/trimmomatic" $PIG_INDEX $NOT_STRAND_SPECIFIC
        echo "INFO: Removing pig genome STAR index from memory" >&2
        $STAR --genomeDir $PIG_INDEX --genomeLoad Remove
elif [ $MONKEY ]; then
	echo "INFO: Loading monkey genome STAR index into memory" >&2
	$STAR --genomeDir $MONKEY_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/trimmomatic" $MONKEY_INDEX $NOT_STRAND_SPECIFIC
	echo "INFO: Removing monkey genome STAR index from memory" >&2
        $STAR --genomeDir $MONKEY_INDEX --genomeLoad Remove
else
	echo "No organism selected. Exiting .... !" >&2
fi


