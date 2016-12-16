#!/bin/bash

#$ -N rna_seq_qc
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/rna_seq_qc.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/rna_seq_qc.out

while getopts ":D:HMPKh" opt; do
        case $opt in
                D)
            rnaseq_dir=$OPTARG
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
        h)
            echo "Usage:    $0 
		-D </absolute/path/to/directory/containing/STAR/output> 
                -H [this sets a flag to use human reference files hg19]
                -M [this sets a flag to use mouse reference files mm9]
		-P [this sets a flag to use pig reference files]
		-K [this sets a flag to use monkey reference files]
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


if [ -z $rnaseq_dir ]; then
	echo "Required Arguments not provided. Try -h to see the documentation"
	exit 1;
fi

if [ -z $HUMAN ] && [ -z $MOUSE ] && [ -z $PIG ] && [ -z $MONKEY ]; then
    echo "ERROR: Either human [-H] or mouse [-M] or pig [-P] or monkey [-K] flags must be provided to the script. Exiting ..... !" >&2
    exit 1
fi

JAVA_EXEC="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/jre1.7.0_80/bin/java"
SAMTOOLS="/usr/bin/samtools"
RNASEQC_EXEC="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/RNA-seq-QC/RNA-SeQC_v1.1.8.jar"
BWA_EXEC="/ifs/rcgroups/mbcf/umv/software/bwa/bwa"

root_dir=$(dirname $rnaseq_dir)
rna_seq_qc_dir="${root_dir}/rnaseq_qc"
mkdir -p $rna_seq_qc_dir
map_file=${rna_seq_qc_dir}/rnaseq_qc_map.txt

if [ -e $map_file ]; then
	rm $map_file;
fi
touch $map_file

echo -e "Sample ID\tBam File\tNotes" >>${map_file}

for sorted_bam_file in $(find $rnaseq_dir -type f -name "*sorted.bam"); do
	id=$(basename $sorted_bam_file | sed -e "s/.sorted.bam//")
	echo -e "${id}\t${sorted_bam_file}\t${id}" >>${map_file}
done

echo "INFO: Processing RNA_seq_QC step" >&2

if [  $MOUSE ]; then
	echo "INFO: Processing mouse reference data" >&2
#	rRNA_ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/mouse_all_rRNA_sequences.fa"
#	ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
#	gtf_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
	ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome_with_ERCC.fa"
	gtf_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Genes/genes_with_ERCC.gtf"
elif [ $HUMAN ]; then
	echo "INFO: Processing human reference data" >&2
	rRNA_ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/RNA-seq-QC/resources/human_rRNA/human_all_rRNA.fasta"
	ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
	gtf_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
	rRNA_interval_list_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/RNA-seq-QC/resources/hg19.rRNA.interval.list"
elif [ $PIG ]; then
	echo "INFO: Processing pig reference data" >&2
	rRNA_ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/RNA-seq-QC/resources/human_rRNA/human_all_rRNA.fasta"
        ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/pig/UCSC/Sus_scrofa/UCSC/susScr3/Sequence/WholeGenomeFasta/genome.fa"
        gtf_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/pig/UCSC/Sus_scrofa/UCSC/susScr3/Annotation/Genes/genes.gtf"
        rRNA_interval_list_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/RNA-seq-QC/resources/hg19.rRNA.interval.list"
elif [ $MONKEY ]; then
	echo "INFO: Processing monkey reference data" >&2
        rRNA_ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/RNA-seq-QC/resources/human_rRNA/human_all_rRNA.fasta"
        ref_fasta="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/monkey/crab_eating/ENSEMBL/Sequence/Macaca_fascicularis_whole_genome.fa"
        gtf_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/monkey/crab_eating/ENSEMBL/Annotation/Macaca_fascicularis.MacFas_5.0.pre.gtf"
        rRNA_interval_list_file="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/RNA-seq-QC/resources/hg19.rRNA.interval.list"
else
	echo "This should never be executed. Some thing wrong ... Exiting!"
	exit 1
fi

$JAVA_EXEC -jar $RNASEQC_EXEC \
	-o ${rna_seq_qc_dir} \
	-r $ref_fasta \
	-s $map_file \
	-t $gtf_file \
	-bwa $BWA_EXEC 
#	-BWArRNA $rRNA_ref_fasta 

