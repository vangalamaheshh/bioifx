#!/bin/bash

# vim: syntax=sh tabstop=4 expandtab

#$ -cwd
#$ -N NS500_preprocess
#$ -e /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/NS500_preprocess.log
#$ -o /zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS/NS500_preprocess.out

function pre_process {
	for dir in $@; do
		if [ -e "${SOURCE_DIR}/${dir}/$RUN_COMPLETION_TOKEN1" ] && [ -e "${SOURCE_DIR}/${dir}/$RUN_COMPLETION_TOKEN2" ]; then
			ILAB_ID=$(grep -iA 2 "^\[Data\]" "${SOURCE_DIR}/${dir}/$RUN_COMPLETION_TOKEN2" | tail -1 | cut -f 1 -d "," | sed -e "s/.*_//")
		fi
		if [ -e "${DEST_DIR}/${dir}" ]; then
			echo "MSG: ${DEST_DIR}/${dir} exists .. so skipping analysis for ${dir}" >&2
			continue
		elif [ -e "${DEST_DIR}/${dir}_${ILAB_ID}" ]; then
			echo "MSG: ${DEST_DIR}/${dir}_${ILAB_ID} exists .. so skipping analysis for ${dir}" >&2
			continue
		else 
			if [ -e "${SOURCE_DIR}/${dir}/$RUN_COMPLETION_TOKEN1" ] && [ -e "${SOURCE_DIR}/${dir}/$RUN_COMPLETION_TOKEN2" ]; then
				echo "MSG: CRON_INFO: Cron job started at $(date)" >&2
				echo "MSG: Processing directory: $dir" >&2
				echo "---------------------------------" >&2
			
				if [ -e "${DEST_DIR}/${dir}/bcl2fastq" ] || [ -e "${DEST_DIR}/${dir}_${ILAB_ID}/bcl2fastq" ]; then
					echo "MSG: bcl2fastq dir exists .. so, skipping bcl2fastq step" >&2
					continue
				else
					mkdir -p ${DEST_DIR}/${dir}_${ILAB_ID}/bcl2fastq

					echo "MSG: Executing '$BASE_MASK_SCRIPT' with run_type=demultiplex_run" >&2
                                        BASE_MASK_INFO=$(/usr/bin/perl $BASE_MASK_SCRIPT --run_info_file="${SOURCE_DIR}/${dir}/RunInfo.xml" --run_type="demultiplex_run")
                                        echo "MSG: Executing Bcl2Fastq run" >&2

					#Now run bcl2fastq script
					command_string="$BCL2FASTQ                                                  \
					--create-fastq-for-index-reads                              \
					--runfolder-dir ${SOURCE_DIR}/${dir}                        \
					--minimum-trimmed-read-length 25                            \
					--no-lane-splitting                                         \
					--input-dir ${SOURCE_DIR}/${dir}/Data/Intensities/BaseCalls \
					--output-dir ${DEST_DIR}/${dir}_${ILAB_ID}/bcl2fastq                   \
					--sample-sheet ${SOURCE_DIR}/${dir}/$RUN_COMPLETION_TOKEN2  \
					--use-bases-mask "${BASE_MASK_INFO}"                        \
					--loading-threads 4                                         \
					--demultiplexing-threads 4                                  \
					--processing-threads 4                                      \
					--writing-threads 4"

					$command_string --barcode-mismatches 1

					if [ $? -ne 0 ]; then
						$command_string --barcode-mismatches 0
					fi

					# remove all index files except for undetermined ones
                                        find ${DEST_DIR}/${dir}_${ILAB_ID}/bcl2fastq -type f -name "[0-9]*_I1_*.fastq.gz" -exec rm {} \;
                                        find ${DEST_DIR}/${dir}_${ILAB_ID}/bcl2fastq -type f -name "[0-9]*_I2_*.fastq.gz" -exec rm {} \;

				fi
			
				if [ -e "${DEST_DIR}/${dir}/concat_per_sample_fastq" ] || [ -e "${DEST_DIR}/${dir}_${ILAB_ID}/concat_per_sample_fastq" ]; then
					echo "MSG: concat_per_sample_fastq dir exists ... so, skipping concat step" >&2
					continue
				else
					mkdir -p ${DEST_DIR}/${dir}_${ILAB_ID}/concat_per_sample_fastq
					ln -s ${DEST_DIR}/${dir}_${ILAB_ID}/concat_per_sample_fastq ${DEST_DIR}/${dir}_${ILAB_ID}/data
					cd ${DEST_DIR}/${dir}_${ILAB_ID}/concat_per_sample_fastq/
					for file in $(find ../bcl2fastq -type f -name "*R1*fastq.gz" -o -name "*R2*fastq.gz"); do
						ln -s $file ${DEST_DIR}/${dir}_${ILAB_ID}/concat_per_sample_fastq/
					done
					date_token=$(echo $dir | grep -Po "^\d+")
                                        rest_nextseq_token=$(echo $dir | grep -Po "_.+?_")
					transfer_dir="${date_token}_${ILAB_ID}${rest_nextseq_token}fastq"
					mkdir -p "$TRANSFER_TO/${transfer_dir}"	
					for file in $(find "${DEST_DIR}/${dir}_${ILAB_ID}/concat_per_sample_fastq" ! -name "Undetermined*"); do
						 cp -f $file $TRANSFER_TO/${transfer_dir}/
					done

					#send mail
					cat > ${MAIL_DIR}/${dir}.NS.mail <<_EOF
NextSeq Run: '${dir}' Status: Completed
Run Type: NextSeq
Run_Name: ${dir}
Completion_time: $(date)
Concat_fastq_files_size:
$(du -hLc ${DEST_DIR}/${dir}_${ILAB_ID}/concat_per_sample_fastq/*.fastq.gz | gawk 'BEGIN{ FS="/"; } { print $1, $NF; }')
_EOF
					#end of send mail
				fi


				if [ -e "${DEST_DIR}/${dir}/fastqc/before_filtering" ] || [ -e "${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/before_filtering" ]; then
					echo "MSG: fastqc/before_filtering dir exists ... so, skipping fastqc step" >&2
					continue
				else
					mkdir -p ${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/before_filtering

					#Now run fastqc step
					FASTQC_INDIR=${DEST_DIR}/${dir}_${ILAB_ID}/bcl2fastq
					#THREADS=$(ls $FASTQC_INDIR/*.fastq.gz | wc -l)
					THREADS=4
					$FASTQC --outdir ${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/before_filtering --threads $THREADS $FASTQC_INDIR/*R*fastq.gz $FASTQC_INDIR/Undetermined*.fastq.gz >&2
				    
				    fastqc_file_list=""
                    for fastqc_file in $(find ${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/before_filtering -type f -name "*_R*_*.zip"); do
                        fastqc_file_list="${fastqc_file_list} -f ${fastqc_file} "
                    done

                    #Generate fastqc stats matrix
                    command="$FASTQC_STATS ${fastqc_file_list}"
                    echo "COMMAND: $command" >&2
                    /usr/bin/perl $command 1>"${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/fastqc_matrix_bf.csv"
                    /zfs/cores/mbcf/mbcf-storage/devel/umv/software/miniconda3/bin/Rscript ${FASTQC_PLOT} "${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/fastqc_matrix_bf.csv" "${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/fastqc_matrix_bf.png" "${DEST_DIR}/${dir}_${ILAB_ID}/fastqc/fastqc_summary_bf.png"

                fi

				echo "MSG: CRON_INFO: Cron job ended at $(date)" >&2
			fi
		fi
	done
}


SOURCE_DIR="/zfs/cores/mbcf/mbcf-storage/NS500/rawdata"
DEST_DIR="/zfs/cores/mbcf/mbcf-storage/NS500/analysis"

TRANSFER_TO="/zfs/cores/mbcf/mbcf-storage/TempTransfer/NS500"

MYHOSTNAME=$(hostname)

#if [[ ! $MYHOSTNAME =~ ^n[0-9]+$ ]]; then 
#	BCL2FASTQ="/ifs/rcgroups/mbcf/umv/software/bcl2fastq/bcl2fastq-build/bin/bcl2fastq"
#	FASTQC="/ifs/rcgroups/mbcf/umv/software/FastQC/fastqc"
#else
	BCL2FASTQ="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/illumina_software/bcl2fastq_v2_17/bin/bcl2fastq"
	FASTQC="/zfs/cores/mbcf/mbcf-storage/devel/umv/software/FastQC/fastqc"
#fi

BASE_MASK_SCRIPT="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/file_types/xml/get_base_mask.pl"

CONCAT_LANE_2_SAMPLE_FASTQ="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/file_types/fastq/concat_lane_2_sample_fastq.pl"
DIR_DIFF_SCRIPT="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/file_types/directories/dir_diff.pl"
FASTQC_STATS="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/stats/fastqc_info.pl"
FASTQC_PLOT="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/stats/fastqc_info.R"

RUN_COMPLETION_TOKEN1="RTAComplete.txt"
RUN_COMPLETION_TOKEN2="SampleSheet.csv"

MAIL_DIR="/zfs/cores/mbcf/mbcf-storage/devel/umv/LOGS"

if [ ! -z "$@" ]; then
	pre_process "$@";
else
	my_dirs=$(/usr/bin/perl $DIR_DIFF_SCRIPT --source_dir=$SOURCE_DIR --dest_dir=$DEST_DIR); 
	pre_process "$my_dirs";
fi

