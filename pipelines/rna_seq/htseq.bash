#!/bin/bash

#$ -N htseq
#$ -pe pvm 10
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/htseq.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/htseq.out

htseq_exec="/apps/python-2.7.9/bin/htseq-count"
htseq_out_dir="/zfs/cores/mbcf/mbcf-storage/NS500/analysis/150830_NS500305_0072_AHJ2J3BGXX/htseq"
MOUSE_GTF_FILE="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf"
mkdir -p $htseq_out_dir

for file in $(find /zfs/cores/mbcf/mbcf-storage/NS500/analysis/150830_NS500305_0072_AHJ2J3BGXX/STAR_with_ERCC_RG_added/ -type f -name "*sorted.bam"); do 
	file_base=$(basename $file | sed -e "s/.sorted.bam/.htseq/"); 
	echo "Processing $file_base" >&2; 
	${htseq_exec} --format=bam --order=pos --stranded=reverse --type=exon --idattr=gene_id --mode=union $file $MOUSE_GTF_FILE 1>${htseq_out_dir}/${file_base}.counts & 
done

wait
echo "All done" >&2
