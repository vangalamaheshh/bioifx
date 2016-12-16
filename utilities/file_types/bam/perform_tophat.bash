#!/bin/bash

#$ -cwd
#$ -N NadjaTophatRun
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/scripts/tophat.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/scripts/tophat.out

human_dir=/zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/hg19
mouse_dir=/zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/mm9

#for file in $(find $human_dir/input -type l); do
#	echo "Processing $file" >&2
#	out_file=$(basename $file)
#	mkdir /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/hg19/output/tophat/$out_file
#	tophat --output-dir /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/hg19/output/tophat/$out_file \
#		--library-type fr-unstranded \
#		--num-threads 12 \
#		--GTF /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf \
#		/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome \
#		$file/*R1*.fastq.gz \
#		$file/*R2*.fastq.gz
#done

for file in $(find $mouse_dir/input -type l); do
        echo "Processing $file" >&2
        out_file=$(basename $file)
	mkdir -p /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/mm9/output/tophat/$out_file
        tophat --output-dir /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/mm9/output/tophat/$out_file \
               --library-type fr-unstranded \
               --num-threads 12 \
               --GTF /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf \
               /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome \
               $file/*R1*.fastq.gz \
               $file/*R2*.fastq.gz
done

exit $?

