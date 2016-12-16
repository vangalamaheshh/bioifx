#!/bin/bash

#$ -N NadjaSTAR
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/scripts/STAR.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/scripts/STAR.out


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
	mouse_dir=/zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/mm9
	for file in $(find $mouse_dir/input -type l); do
		echo "Processing $file" >&2
	        out_file=$(basename $file)
		mkdir -p /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/mm9/output/STAR/$out_file
		$STAR --runMode alignReads \
			--runThreadN 12 \
			--genomeDir $STAR_INDEX_DIR \
			--readFilesIn $file/*R1*.fastq.gz $file/*R2*.fastq.gz \
			--readFilesCommand zcat \
			--outFileNamePrefix /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/mm9/output/STAR/$out_file/$out_file \
			--outSAMmode Full \
			--outSAMattributes All
	done

}

function remove_STAR_index {

	# Remove loaded STAR index from memory
	$STAR --genomeDir $STAR_INDEX_DIR --genomeLoad Remove

}

#STAR="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/STAR/STAR_2.3.0e.Linux_x86_64_static/STAR"
STAR="/usr/local/bin/STAR"
STAR_INDEX_DIR="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/STARIndex"
GENOME_PATH="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
GTF_PATH="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"

generate_STAR_index
load_STAR_index
map_reads
remove_STAR_index

exit $?

