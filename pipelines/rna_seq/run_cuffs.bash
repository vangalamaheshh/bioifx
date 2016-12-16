#!/bin/bash
#$ -N RUN_CUFFS
#$ -pe pvm 8
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_cuffs.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/run_cuffs.out


CUFFLINKS="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/cufflinks-2.2.1.Linux_x86_64/cufflinks"
Monkey_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/monkey/crab_eating/Annotation/CE.gtf"
HUMAN_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes_with_ERCC.gtf"
MOUSE_GTF="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf"

for sorted_bam in $(find /zfs/cores/mbcf/mbcf-storage/NS500/analysis/151030_NS500233_0255_AHL7KVBGXX_AT2356/STAR_hg19/ -type f -name "*sorted.bam"); do
	out_prefix=$(basename $sorted_bam | sed -e "s/.sorted.bam//")
	out_dir="/zfs/cores/mbcf/mbcf-storage/NS500/analysis/151030_NS500233_0255_AHL7KVBGXX_AT2356/cufflinks_hg19/${out_prefix}"
	mkdir -p ${out_dir}
	$CUFFLINKS -o ${out_dir} -p 8 -G $HUMAN_GTF $sorted_bam
	for cuff_out_file in $(find ${out_dir} -mindepth 1 -maxdepth 1 -type f ! -name "\.*"); do
		rename_file="${out_prefix}.$(basename $cuff_out_file)"
		mv $cuff_out_file "${out_dir}/${rename_file}"
	done
done

echo "Done" >&2

