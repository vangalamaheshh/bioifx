#!/bin/bash
# vim: syntax=sh tabstop=4 expandtab

#---------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Sep, 12, 2016
#----------------------------------

set -euo pipefail

analysis_dir=$1

fastqc_file_list=""
for fastqc_file in $(find ${analysis_dir}/fastqc/before_filtering -type f -name "*_R*_*.zip"); do
    fastqc_file_list="${fastqc_file_list} -f ${fastqc_file} "
done

#Generate fastqc stats matrix
FASTQC_STATS="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/stats/fastqc_info.pl"
FASTQC_PLOT="/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/utilities/stats/fastqc_info.R"
command="$FASTQC_STATS ${fastqc_file_list}"
echo "COMMAND: $command" >&2
/usr/bin/perl $command 1>"${analysis_dir}/fastqc/fastqc_matrix_bf.csv"
/zfs/cores/mbcf/mbcf-storage/devel/umv/software/miniconda3/bin/Rscript ${FASTQC_PLOT} "${analysis_dir}/fastqc/fastqc_matrix_bf.csv" "${analysis_dir}/fastqc/fastqc_matrix_bf.png" \
 "${analysis_dir}/fastqc/fastqc_summary_bf.png"

exit $?
