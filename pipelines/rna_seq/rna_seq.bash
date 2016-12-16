#!/bin/bash

#$ -N RNA_SEQ
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/rna_seq.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/rna_seq.out


function generate_STAR_index {
	
	$STAR --runMode genomeGenerate --genomeDir $STAR_INDEX_DIR --genomeFastaFiles $GENOME_PATH \
		--runThreadN 8 --sjdbGTFfile $GTF_PATH --sjdbOverhang 99

}

function load_STAR_index {
	
	# Load star index into memory and keep it for future runs
	$STAR --genomeDir $STAR_INDEX_DIR --genomeLoad LoadAndExit

}

function map_reads {

	# Run mapping for HIV
	$STAR --runMode alignReads \
		--runThreadN 8 \
		--genomeDir $STAR_INDEX_DIR \
		--readFilesIn /zfs/cores/mbcf/mbcf-storage/devel/umv/test/test_rna_seq/input/081114_HIV_Engelman_TGCAG_L001_R1_001.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix /zfs/cores/mbcf/mbcf-storage/devel/umv/test/test_rna_seq/output/HIV/HIV \
		--outSAMmode Full \
		--outSAMattributes All

	# Run mapping for MxB
	$STAR --runMode alignReads \
		--runThreadN 8 \
		--genomeDir $STAR_INDEX_DIR \
		--readFilesIn /zfs/cores/mbcf/mbcf-storage/devel/umv/test/test_rna_seq/input/081114_MxB_Engelman_ACGTC_L001_R1_001.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix /zfs/cores/mbcf/mbcf-storage/devel/umv/test/test_rna_seq/output/MxB/MxB \
		--outSAMmode Full \
		--outSAMattributes All

}

function remove_STAR_index {

	# Remove loaded STAR index from memory
	$STAR --genomeDir $STAR_INDEX_DIR --genomeLoad Remove

}

STAR="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR_2.3.0e.Linux_x86_64_static/STAR"
STAR_INDEX_DIR="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human_hg19/STAR_index"
GENOME_PATH="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human_hg19/genome.fa"
GTF_PATH="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human_hg19/genes.gtf"

#generate_STAR_index
load_STAR_index
map_reads
remove_STAR_index

exit $?

