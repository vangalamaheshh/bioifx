#!/bin/bash

# Takes input_directory path, output_directory path, output_prefix

INPUT_DIR=$1
OUTPUT_DIR=$2

INPUT_DIR=$(echo $INPUT_DIR | sed -e "s/\/$//")
OUTPUT_DIR=$(echo $OUTPUT_DIR | sed -e "s/\/$//")

FASTQC="/media/mvangala/data/software/bioinfo/QC/FastQC/fastqc"
FASTQC_RAW_DIR="fastqc/raw_reads"
FASTQC_FILTERED_DIR="fastqc/filtered_reads"

TRIMMOMATIC="/media/mvangala/data/software/bioinfo/QC/Trimmomatic-0.32/trimmomatic-0.32.jar"
ILLUMINA_ADAPTER_FILE="/media/mvangala/data/software/bioinfo/QC/Trimmomatic-0.32/adapters/TruSeq3-SE.fa"
FILTERED_READ_DIR="filtered_reads"

FASTQ_2_FASTA="/media/mvangala/data/software/bioinfo/format/fastx_toolkit-0.0.14/src/fastq_to_fasta/fastq_to_fasta"
FASTQ_2_FASTA_OUT_DIR="fastq_2_fasta"

UCLUST="/media/mvangala/data/software/bioinfo/alignment/uclust/uclustq1.2.22_i86linux64"
UCLUST_OUT_DIR="uclust"

FETCH_UCLUST_UNIQ_SEQS="/media/mvangala/mbcf_storage_aberdeen/Engelman/scripts/fetch_uclust_uniq_seqs.pl"

PYTHON3="/usr/bin/python3"
FETCH_SEQS_WITH_LINKER="/media/mvangala/mbcf_storage_aberdeen/Engelman/scripts/Hiv_seq_match.py"
LINKER_OUT_DIR="LTR"

BOWTIE="/media/mvangala/data/software/bioinfo/alignment/bowtie/bowtie"
GENOME_INDEX="/media/mvangala/mbcf_storage_aberdeen/ref_files/human/Homo_sapiens/Ensembl/GRCh38/Sequences/BowtieIndex/GRCh_genome"
BOWTIE_OUT_DIR="bowtie"

HTSEQ="/media/mvangala/mbcf_storage_aberdeen/Engelman/scripts/HTSeq_gene_stats.py"
GTF_FILE="/media/mvangala/mbcf_storage_aberdeen/ref_files/human/Homo_sapiens/Ensembl/GRCh38/Annotation/Homo_sapiens.GRCh38.76.gtf"
HTSEQ_OUT_DIR="htseq"
PYTHON2="/usr/bin/python"

function run_fastqc {
	local fastqc=$1
	local in_dir=$2
	local out_dir=$3
	
	mkdir -p $out_dir
	$fastqc --outdir $out_dir --threads $(ls $in_dir/*fastq.gz | wc -l) ${in_dir}/*fastq.gz
}

function run_trimmomatic {
	local trimmomatic=$1
	local in_dir=$2
	local out_dir=$3
    local illumina_adapter_file=$4

	mkdir -p $out_dir
	for file in ${in_dir}/*fastq.gz ; do
		local out_file=$(basename $file | sed -e "s/.fastq.gz$/.filtered.fastq.gz/")
		java -jar $trimmomatic SE -threads 8 -phred33 $file "${out_dir}/${out_file}" ILLUMINACLIP:${ILLUMINA_ADAPTER_FILE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done
}

function run_fastq_2_fasta {
	local fastq_2_fasta=$1
	local in_dir=$2
	local out_dir=$3
	
	mkdir -p $out_dir
	for file in ${in_dir}/*fastq.gz ; do
		local out_file=$(basename $file | sed -e "s/.fastq.gz$/.fasta/")
		$fastq_2_fasta -i <(zcat $file) -o "${out_dir}/${out_file}" &
	done
	wait
}

function run_uclust {
	local uclust=$1
	local in_dir=$2
	local out_dir=$3
	
	mkdir -p $out_dir
	for file in ${in_dir}/*.fasta ; do
		local sorted_file=$(basename $file | sed -e "s/.fasta$/.sorted.fasta/")
		local uclust_out_file=$(basename $file | sed -e "s/.filtered.fasta$/.uclust.out/")
		$uclust --mergesort $file --output "${out_dir}/${sorted_file}"
		$uclust --input "${out_dir}/${sorted_file}" --uc "${out_dir}/${uclust_out_file}" --id 0.95
	done
}

function run_fetch_uclust_uniq_seqs {
	local fetch_uclust_uniq_seqs=$1
	local uclust_dir=$2
	perl $fetch_uclust_uniq_seqs --uclust_dir $uclust_dir
}

function run_fetch_seqs_with_linker {
	local python3=$1
	local linker_script=$2
	local in_dir=$3
	local out_dir=$4
	
	for file in ${in_dir}/*.uclust.uniq.fasta ; do
		local token=$(basename $file | sed -e "s/.uclust.uniq.fasta//")
		mkdir -p "${out_dir}/${token}"
		$python3 $linker_script --file $file --output_token "${out_dir}/${token}/${token}"
	done
}

function run_bowtie {
    local bowtie=$1
    local genome_index=$2
    local in_dir=$3
    local out_dir=$4
    
    mkdir -p ${out_dir}
    for linker_dir in $(find ${in_dir} -maxdepth 1 -mindepth 1 -type d); do
        local token=$(basename ${linker_dir})
        $bowtie -l 25 -t --sam --threads 8 -f -m 1 $genome_index <(zcat ${linker_dir}/*.LTR_and_LINKER.fasta.gz ${linker_dir}/*.LTR_but_no_LINKER.fasta.gz) 1>"${out_dir}/${token}.bowtie_aligned.sam"
    done
}

function run_htseq {
    local python2=$1
    local htseq=$2
    local gtf_file=$3
    local in_dir=$4
    local out_dir=$5

    mkdir -p ${out_dir}
    for sam_file in ${in_dir}/*.bowtie_aligned.sam ; do
        local token=$(basename $sam_file | sed -e "s/.bowtie_aligned.sam//")
        $python2 $htseq $sam_file $gtf_file 1>"${out_dir}/${token}.htseq.out"
    done
}

function log_info {
    local token=$1
    echo ""
}

run_fastqc $FASTQC "$INPUT_DIR" "${OUTPUT_DIR}/${FASTQC_RAW_DIR}"
run_trimmomatic $TRIMMOMATIC "$INPUT_DIR" "${OUTPUT_DIR}/${FILTERED_READ_DIR}" $ILLUMINA_ADAPTER_FILE
run_fastqc $FASTQC "${OUTPUT_DIR}/${FILTERED_READ_DIR}" "${OUTPUT_DIR}/${FASTQC_FILTERED_DIR}"
run_fastq_2_fasta $FASTQ_2_FASTA "${OUTPUT_DIR}/${FILTERED_READ_DIR}" "${OUTPUT_DIR}/${FASTQ_2_FASTA_OUT_DIR}"
run_uclust $UCLUST "${OUTPUT_DIR}/${FASTQ_2_FASTA_OUT_DIR}" "${OUTPUT_DIR}/${UCLUST_OUT_DIR}"
run_fetch_uclust_uniq_seqs $FETCH_UCLUST_UNIQ_SEQS "${OUTPUT_DIR}/${UCLUST_OUT_DIR}"
run_fetch_seqs_with_linker $PYTHON3 $FETCH_SEQS_WITH_LINKER "${OUTPUT_DIR}/${UCLUST_OUT_DIR}" "${OUTPUT_DIR}/${LINKER_OUT_DIR}"
run_bowtie $BOWTIE $GENOME_INDEX "${OUTPUT_DIR}/${LINKER_OUT_DIR}" "${OUTPUT_DIR}/${BOWTIE_OUT_DIR}"
run_htseq $PYTHON2 $HTSEQ $GTF_FILE "${OUTPUT_DIR}/${BOWTIE_OUT_DIR}" "${OUTPUT_DIR}/${HTSEQ_OUT_DIR}"
exit $?;
