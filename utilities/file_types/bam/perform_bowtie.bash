#!/bin/bash

#$ -cwd
#$ -N NadjaBowiteRun
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/scripts/bowtie.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/scripts/bowtie.out

human_dir=/zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/hg19
mouse_dir=/zfs/cores/mbcf/mbcf-storage/devel/umv/nadja_bowtie/mm9

#for file in $(find $human_dir/input -type l); do
#	echo "Processing $file" >&2
#	out_file=$(basename $file)
#	/zfs/cores/mbcf/mbcf-storage/devel/umv/software/bowtie/bowtie-1.0.1/bowtie -t --threads 12 /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome -1 <(zcat $file/*R1*.fastq.gz) -2 <(zcat $file/*R2*.fastq.gz) 1>$human_dir/output/$out_file.sam.out 2>$human_dir/output/$out_file.log
#done

for file in $(find $mouse_dir/input -type l); do
        echo "Processing $file" >&2
        out_file=$(basename $file)
        /zfs/cores/mbcf/mbcf-storage/devel/umv/software/bowtie/bowtie-1.0.1/bowtie -t --threads 12 /zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome -1 <(zcat $file/*R1*.fastq.gz) -2 <(zcat $file/*R2*.fastq.gz) 1>$mouse_dir/output/$out_file.sam.out 2>$mouse_dir/output/$out_file.log
done

exit $?

