#!/bin/bash

#$ -N RSeQC
#$ -pe pvm 4
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/rseqc.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/rseqc.out

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

RSEQC_GB_EXEC="/apps/python-2.7.9/bin/geneBody_coverage.py"
RSEQC_RD_EXEC="/apps/python-2.7.9/bin/read_distribution.py"

root_dir=$(dirname $rnaseq_dir)
rseqc_dir="${root_dir}/rseqc"
mkdir -p $rseqc_dir
mkdir -p ${rseqc_dir}/gene_body
mkdir -p ${rseqc_dir}/read_dist

sorted_bam_files=""
for sorted_bam_file in $(find $rnaseq_dir -type f -name "*sorted.bam"); do
	if [ -z $sorted_bam_files ]; then
		sorted_bam_files="${sorted_bam_file}"
		ilab_id=$(basename $sorted_bam_file | sed -e "s/.sorted.bam//" | perl -e 'my $ilab_id = <STDIN>; chomp $ilab_id; my @ilab = split("-", $ilab_id ); print join("-", @ilab[2 .. scalar @ilab - 1]);')
	else
		sorted_bam_files="${sorted_bam_files},${sorted_bam_file}"
	fi
done

echo "INFO: Processing Gene Body Coverage step" >&2

if [  $MOUSE ]; then
	echo "INFO: Processing mouse reference data" >&2
	ref_gene_model="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/mouse/Mus_musculus/UCSC/mm9/Annotation/Genes/mm9_genes.bed"
elif [ $HUMAN ]; then
	echo "INFO: Processing human reference data" >&2
	ref_gene_model="/zfs/cores/mbcf/mbcf-storage/devel/umv/ref_files/human/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19_genes.bed"
elif [ $PIG ]; then
	echo "INFO: Processing pig reference data" >&2
elif [ $MONKEY ]; then
	echo "INFO: Processing monkey reference data" >&2
else
	echo "This should never be executed. Some thing wrong ... Exiting!"
	exit 1
fi

 $RSEQC_GB_EXEC -i ${sorted_bam_files} \
	-r $ref_gene_model \
	-f png \
	-o ${rseqc_dir}/gene_body/${ilab_id}


echo "INFO: Processing Read Distribution step" >&2

i=1; j=4;
for file in $(echo $sorted_bam_files | tr "," " "); do
	out_file=$(basename $file | sed -e "s/.sorted.bam//");
	if [ $i -lt $j ]; then
		i=$[ $i + 1 ]
		echo "Processing $file" >&2
		$RSEQC_RD_EXEC -i $file -r $ref_gene_model 1>${rseqc_dir}/read_dist/${out_file}.txt &		 
	else
		echo "Processing $file" >&2
		$RSEQC_RD_EXEC -i $file -r $ref_gene_model 1>${rseqc_dir}/read_dist/${out_file}.txt &
		i=1
		wait
	fi
done

wait



