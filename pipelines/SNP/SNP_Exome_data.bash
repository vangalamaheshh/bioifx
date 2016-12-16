#!/bin/bash

ANALYSIS_DIR="/media/mvangala/EXOME_DATA_150729_NS500144_0302_AHGLHYBGXX/SNP"
REF_FASTA="/media/mvangala/mbcf_storage_aberdeen/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
STAR_GENOME_INDEX_DIR="/media/mvangala/mbcf_storage_aberdeen/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/STARIndex"

function trim_reads {
	local left_file=$1
	local right_file=$2
	local trim_dir=$3
	mkdir -p $trim_dir
	filtered_left_paired=$(basename $left_file | sed -e "s/.fastq/.trimmomatic.left_paired.fastq/")
	filtered_left_unpaired==$(basename $left_file | sed -e "s/.fastq/.trimmomatic.left_unpaired.fastq/")
	filtered_right_paired=$(basename $right_file | sed -e "s/.fastq/.trimmomatic.right_paired.fastq/")
	filtered_right_unpaired=$(basename $right_file | sed -e "s/.fastq/.trimmomatic.right_unpaired.fastq/")
	java -jar /media/mvangala/data/software/bioinfo/QC/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 \
		$left_file $right_file \
		"${trim_dir}/${filtered_left_paired}" "${trim_dir}/${filtered_left_unpaired}" \
       	"${trim_dir}/${filtered_right_paired}" "${trim_dir}/${filtered_right_unpaired}" \
      	ILLUMINACLIP:/media/mvangala/data/software/bioinfo/QC/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 \
       	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
}

function align_reads {
	local trim_dir=$1
	local align_dir=$2
	local index_dir=$3
	mkdir -p $align_dir
	left_file=$(find $trim_dir -type f -name "*.trimmomatic.left_paired.fastq.gz")
	right_file=$(find $trim_dir -type f -name "*.trimmomatic.right_paired.fastq.gz")
	cd $align_dir
	/media/mvangala/data/software/bioinfo/alignment/STAR-RNA/STAR_2.3.0e/STAR --genomeDir $index_dir --readFilesIn $left_file $right_file \
		--runThreadN 6 --readFilesCommand zcat 
}

function generate_index {
	local two_pass_index_dir=$1
	local ref_fasta=$2
	local sjdb_file=$3
	mkdir -p ${two_pass_index_dir}
	/media/mvangala/data/software/bioinfo/alignment/STAR-RNA/STAR_2.3.0e/STAR --runMode genomeGenerate --genomeDir ${two_pass_index_dir} \
		--genomeFastaFiles $ref_fasta \
		--sjdbFileChrStartEnd $sjdb_file \
		--sjdbOverhang 75 --runThreadN 6
}

function process_picard {
	local two_pass_dir=$1
	local picard_dir=$2
	local sample_name=$3
	mkdir -p $picard_dir
	#Add or replace read groups
	java -jar /media/mvangala/data/software/bioinfo/format/picard-tools-1.115/AddOrReplaceReadGroups.jar I="${two_pass_dir}/Aligned.out.sam" \
		O="${picard_dir}/${sample_name}.rg_added_sorted.bam" SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample

	#Mark duplicates
	java -jar /media/mvangala/data/software/bioinfo/format/picard-tools-1.115/MarkDuplicates.jar I="${picard_dir}/${sample_name}.rg_added_sorted.bam" \
		O="${picard_dir}/${sample_name}.dedupped.bam" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
		M="${picard_dir}/${sample_name}.output.metrics"
}

function split_n_trim_reads {
	local ref_fasta=$1
	local picard_dir=$2
	local sample_name=$3
	local split_n_trim_dir=$4
	mkdir -p $split_n_trim_dir
	java -jar /media/mvangala/data/software/bioinfo/frameworks/GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R $ref_fasta \
		-I "${picard_dir}/${sample_name}.dedupped.bam" -o "${split_n_trim_dir}/${sample_name}.split.bam" \
		-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
}

function perform_variant_calling {
	local ref_fasta=$1
	local split_n_trim_dir=$2
	local sample_name=$3
	local variant_calling_dir=$4
	mkdir -p $variant_calling_dir
	#variant calling
	java -jar /media/mvangala/data/software/bioinfo/frameworks/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref_fasta \
		-I "${split_n_trim_dir}/${sample_name}.split.bam" -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 \
		-o "${variant_calling_dir}/${sample_name}.output.vcf"

	#variant filtering
	java -jar /media/mvangala/data/software/bioinfo/frameworks/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $ref_fasta \
		-V "${variant_calling_dir}/${sample_name}.output.vcf" -window 35 -cluster 3 -o "${variant_calling_dir}/${sample_name}.filtered_output.vcf"
}

for left_file in $(find /media/mvangala/EXOME_DATA_150729_NS500144_0302_AHGLHYBGXX/ -name "[0-9]*R1*fastq.gz" -a ! -iname "*Exome*"); do 
	echo "INFO:Processing sample: $left_file" >&2
	echo "-----------------------------------" >&2
	right_file=$(echo $left_file | sed -e "s/_R1/_R2/")
	sample_base_dir=$(basename $(dirname $(dirname $left_file)))
	sample_name=$(basename $left_file | sed -e "s/_.*//");
#	sample_name=${sample_base_dir}
#	sample_base_dir="${sample_base_dir}/${sample_name}"
#	trim_reads $left_file $right_file "${ANALYSIS_DIR}/${sample_name}/trimmomatic"
#	align_reads "${ANALYSIS_DIR}/${sample_name}/trimmomatic" "${ANALYSIS_DIR}/${sample_name}/1pass" ${STAR_GENOME_INDEX_DIR}
#	generate_index "${ANALYSIS_DIR}/${sample_name}/2pass_index" $REF_FASTA "${ANALYSIS_DIR}/${sample_name}/1pass/SJ.out.tab" 
#	align_reads "${ANALYSIS_DIR}/${sample_name}/trimmomatic" "${ANALYSIS_DIR}/${sample_name}/2pass" "${ANALYSIS_DIR}/${sample_name}/2pass_index"
#	process_picard "${ANALYSIS_DIR}/${sample_name}/2pass" "${ANALYSIS_DIR}/${sample_name}/picard" $sample_name
#	split_n_trim_reads $REF_FASTA  "${ANALYSIS_DIR}/${sample_name}/picard" $sample_name "${ANALYSIS_DIR}/${sample_name}/split_n_trim"
	perform_variant_calling $REF_FASTA "${ANALYSIS_DIR}/${sample_name}/split_n_trim" $sample_name "${ANALYSIS_DIR}/${sample_name}/variant_calling"
done




































