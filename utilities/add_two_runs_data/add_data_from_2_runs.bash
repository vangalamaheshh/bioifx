#!/bin/bash

set -euo pipefail

old_run="160713_NB501431_0106_AHF2JLBGXY_EK3174/data/"
new_run="160716_NS500144_0569_AHGHYNBGXY_EK3174/data/"

out_run="160716_NB501431_0106_CONCAT_AHGHYNBGXY_EK3174/data/"

for file in $(find $old_run -type l -name "20*_R2_*"); do
	token=$(basename $file | grep -oP ".+\-EK3174")
	echo "Processing $token" >&2
	new_run_file=$(find $new_run -type l -name "*${token}*_R2_*")
	out_name=$(basename $file | sed -e "s/20160713/20160713/g")
	cat $file $new_run_file 1>${out_run}/${out_name}
done
