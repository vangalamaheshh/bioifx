#!/bin/bash

set -e 

my_date=$(date +"%m-%d-%y_%H:%M:%S")
my_date="${my_date}_$$"
my_tag_name="comparative_track_${my_date}"
my_wrapper_name="comparative_wrapper_${my_date}"
my_worker_name="comparative_worker_${my_date}"
my_temp_file="/tmp/comparative_track_${my_date}.config"

track_name=$1
acc_ids=$2

umask 000

touch "$my_temp_file"

#...............................................
# This work around here will create an empty tag.
# Also, if the tag exists, it doesn't change those
#+files
#-----------------------------------------------
vp-add-dataset --tag-name ${my_tag_name} --append
#-----------------------------------------------

vp_describe_command="vp-describe-protocols -p ${track_name} -c input.ACC_IDS=${acc_ids} -c input.GENBANK_TAG=${my_tag_name} -c pipeline.PIPELINE_NAME=${my_worker_name}"

#echo $vp_describe_command

$vp_describe_command 1>$my_temp_file

run_pipeline="vp-run-pipeline --print-task-name --pipeline-config $my_temp_file"

#echo $run_pipeline

task_name=$($run_pipeline)

if [ "$?" -ne 0 ]; then
	echo 'vp-run-pipeline failed to run'
	exit 1
fi

echo $task_name

exit $?
