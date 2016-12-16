#!/bin/bash

#$ -N RUN_STAR
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_STAR.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_STAR.out
#$ -cwd
#$ -pe pvm 8
#$ -l mem_free=45G 
#$ -q all.q

function map_reads {
	local input_dir=$1
	local index_dir=$2
	local not_strand_specific=$3
	local no_gene_counts=$4
	local run_trim=$5
	local output_dir=""
	local trim_dir=""
	if [ $run_trim = 1 ]; then
		output_dir="$(dirname ${input_dir})/STAR_Trim"
		trim_dir="$(dirname ${input_dir})/trimmomatic"
	else
		output_dir="$(dirname ${input_dir})/STAR"
	fi
	local strand_command=""
	local gene_counts_command=""

	if [ $not_strand_specific = 1 ]; then
		echo "INFO: Processing not strand specific STAR" >&2
		strand_command="--outSAMstrandField intronMotif"
	else
		echo "INFO: Processing strand specific STAR" >&2
		strand_command="--outFilterIntronMotifs RemoveNoncanonical"
	fi
	
	if [ $no_gene_counts = 1 ]; then
		echo "INFO: Gene counts feature will not be processed" >&2
	else
		gene_counts_command="--quantMode GeneCounts"
	fi

	for left_mate in $(find $input_dir -mindepth 1 -maxdepth 1 ! -name "Undetermined*fastq.gz" | grep ".*_R1_.*.fastq.gz"); do
		right_mate=$(echo $left_mate | sed -e "s/_R1/_R2/")
		output_token=${META_INFO[$(basename ${left_mate})]} 
             	echo "INFO: Processing ${output_token}" >&2
		if [ $run_trim = 1 ]; then	
			echo "INFO: Processing trimmomatic" >&2
			$TRIM_SCRIPT $output_token ${trim_dir}/${output_token} $left_mate $right_mate
			if [ $? -ne 0 ]; then
				echo "Trimmomatic step exited with a non-zero error. Exiting .... !" >&2
				exit 1;
			fi
			left_mate="${trim_dir}/${output_token}/${output_token}.trimmomatic.left_paired.fastq.gz"
			right_mate="${trim_dir}/${output_token}/${output_token}.trimmomatic.right_paired.fastq.gz"
		fi
		mkdir -p ${output_dir}/${output_token}
		if [ ! -f $right_mate ]; then
			local command="$STAR --runMode alignReads --runThreadN 8 --genomeDir $index_dir --readFilesIn $left_mate --readFilesCommand zcat --outFileNamePrefix ${output_dir}/${output_token}/${output_token}. --outSAMmode Full --outSAMattributes All  ${strand_command} --outSAMattrRGline ID:${output_token} PL:illumina LB:${output_token} SM:${output_token} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 45000000000 ${gene_counts_command}"	
			echo "INFO: Executing command (SE): $command" >&2
			$command
		else
			local command="$STAR --runMode alignReads --runThreadN 8 --genomeDir $index_dir --readFilesIn $left_mate $right_mate --readFilesCommand zcat --outFileNamePrefix ${output_dir}/${output_token}/${output_token}. --outSAMmode Full --outSAMattributes All $strand_command --outSAMattrRGline ID:${output_token} PL:illumina LB:${output_token} SM:${output_token} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 45000000000 ${gene_counts_command}"
			echo "INFO: Executing command (PE): $command" >&2
			$command
		fi
		mv "${output_dir}/${output_token}/${output_token}.Aligned.sortedByCoord.out.bam" "${output_dir}/${output_token}/${output_token}.sorted.bam"
		if [ $no_gene_counts = 0 ]; then
			mv "${output_dir}/${output_token}/${output_token}.ReadsPerGene.out.tab" "${output_dir}/${output_token}/${output_token}.counts.tab"
		fi
		$SAMTOOLS index "${output_dir}/${output_token}/${output_token}.sorted.bam"
	done
	
	wait
	mkdir -p ${output_dir}/BAM
	cd ${output_dir}/BAM
	for file in $(find ../../$(basename ${output_dir})/ -type f -name "*sorted.ba*" ); do ln -s $file ${output_dir}/BAM/ ; done	
	STAR_LOG_FILES=""
	STAR_COUNT_FILES=""
	for sample_name in ${SAMPLE_LIST[@]}; do 
		STAR_LOG_FILES="${STAR_LOG_FILES} -l ${output_dir}/${sample_name}/${sample_name}.Log.final.out "; 
		STAR_COUNT_FILES="${STAR_COUNT_FILES} -f ${output_dir}/${sample_name}/${sample_name}.counts.tab ";
	done
	perl $STAR_REPORT $STAR_LOG_FILES 1>${output_dir}/STAR_Align_Report.csv
	Rscript $MAP_STATS ${output_dir}/STAR_Align_Report.csv ${output_dir}/STAR_Align_Report.png
	if [ $not_strand_specific ] && [ $no_gene_counts = 0 ]; then
		perl $STAR_COUNT_MATRIX --column 1 $STAR_COUNT_FILES 1>${output_dir}/STAR_Gene_Counts.csv
	elif [ $no_gene_counts = 0 ]; then
		perl $STAR_COUNT_MATRIX $STAR_COUNT_FILES 1>${output_dir}/STAR_Gene_Counts.csv
	fi
}

STAR="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR"
STAR_COUNT_MATRIX="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/rna_seq/raw_and_fpkm_count_matrix.pl"
STAR_REPORT="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/report_generation/STAR_reports.pl"
MAP_STATS="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/report_generation/map_stats.R"
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

TRIM_SCRIPT="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/QC/trimmomatic/perform_trimmomatic.bash"

SAMTOOLS="/usr/bin/samtools"

NOT_STRAND_SPECIFIC=0
NO_GENE_COUNTS=0
RUN_TRIM=0

while getopts ":D:HMPKSF:I:GTh" opt; do
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
		F)
			META_SHEET=$OPTARG
			;;
		I)
			USER_STAR_INDEX=$OPTARG
			;;
		G)
			NO_GENE_COUNTS=1
			;;
		T)
			RUN_TRIM=1
			;;
		h)
			echo "Usage: 	$0 
				-D <directory>
				-F <meta sheet> 
				-H [this sets a flag to use human STAR index hg19]
				-P [this sets a flag to use pig STAR index]
				-M [this sets a flag to use mouse STAR index mm9]
				-K [this sets a flag to use monkey STAR index]
				-S [if you run STAR as not strand-specific; use this flag]
				-I [if you have STAR index for any other organism, use this flag to provide absolute path to that index]
				-G [if you want to run STAR with out gene counts feature, use this flag]
				-T [if you want to run trimmomatic before running STAR, use this flag]
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

if [ -z $HUMAN ] && [ -z $MOUSE ] && [ -z $PIG ] && [ -z $MONKEY ] && [ -z $USER_STAR_INDEX ]; then
	echo "ERROR: Either human [-H] or mouse [-M] or pig [-P] or monkey [-K] flags or a custom STAR index must be provided to the script. Exiting ..... !" >&2
	exit 1
fi

if [ -z $META_SHEET ]; then
	echo "ERROR: User did not provide a metasheet. Exiting .... !" >&2
	exit 1
fi

declare -A META_INFO
declare -a SAMPLE_LIST

for line in $(sed -ne "2,$ p" $META_SHEET); do
	IFS=',' read -a my_array <<< "$line"
	SAMPLE_LIST=(${SAMPLE_LIST[@]} ${my_array[1]})
	META_INFO[${my_array[0]}]=${my_array[1]}
done

cd #change to home since STAR creates useless files in the cwd

if [[ $RUN_DIR == *MISEQ* ]]; then
	FASTQ_IN_DIR="/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/${RUN_DIR}/concat_per_sample_fastq"
else
	FASTQ_IN_DIR="/zfs/cores/mbcf/mbcf-storage/NS500/analysis/${RUN_DIR}/concat_per_sample_fastq"
fi

if [ ! -d $FASTQ_IN_DIR ]; then
	FASTQ_IN_DIR="/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/${RUN_DIR}/concat_per_sample_fastq"
fi

if [ $HUMAN ]; then
	echo "INFO: Loading human genome STAR index into memory" >&2
	$STAR --genomeDir $HG19_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads $FASTQ_IN_DIR $HG19_INDEX $NOT_STRAND_SPECIFIC $NO_GENE_COUNTS $RUN_TRIM
	echo "INFO: Removing human genome STAR index from memory" >&2
	$STAR --genomeDir $HG19_INDEX --genomeLoad Remove
elif [ $MOUSE ]; then
	echo "INFO: Loading mouse genome STAR index into memory" >&2
	$STAR --genomeDir $MM9_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads $FASTQ_IN_DIR $MM9_INDEX $NOT_STRAND_SPECIFIC $NO_GENE_COUNTS $RUN_TRIM
	echo "INFO: Removing mouse genome STAR index from memory" >&2
	$STAR --genomeDir $MM9_INDEX --genomeLoad Remove
elif [ $PIG ]; then
	echo "INFO: Loading pig genome STAR index into memory" >&2
        $STAR --genomeDir $PIG_INDEX --genomeLoad LoadAndKeep
        echo "INFO: Mapping reads" >&2
        map_reads $FASTQ_IN_DIR $PIG_INDEX $NOT_STRAND_SPECIFIC $NO_GENE_COUNTS $RUN_TRIM
        echo "INFO: Removing pig genome STAR index from memory" >&2
        $STAR --genomeDir $PIG_INDEX --genomeLoad Remove
elif [ $MONKEY ]; then
	echo "INFO: Loading monkey genome STAR index into memory" >&2
	$STAR --genomeDir $MONKEY_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads $FASTQ_IN_DIR $MONKEY_INDEX $NOT_STRAND_SPECIFIC $NO_GENE_COUNTS $RUN_TRIM
	echo "INFO: Removing monkey genome STAR index from memory" >&2
        $STAR --genomeDir $MONKEY_INDEX --genomeLoad Remove
elif [ $USER_STAR_INDEX ]; then
	echo "INFO: Loading user provided STAR index into memory" >&2
	$STAR --genomeDir $USER_STAR_INDEX --genomeLoad LoadAndKeep
	echo "INFO: Mapping reads" >&2
	map_reads $FASTQ_IN_DIR $USER_STAR_INDEX $NOT_STRAND_SPECIFIC $NO_GENE_COUNTS $RUN_TRIM
	echo "INFO: Removing user specificed genome STAR index from memory" >&2
        $STAR --genomeDir $USER_STAR_INDEX --genomeLoad Remove
else
	echo "No organism selected. Exiting .... !" >&2
fi


