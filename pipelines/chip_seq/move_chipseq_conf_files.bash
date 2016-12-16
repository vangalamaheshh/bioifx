#!/bin/bash

read -a conf_file_list <<< $(find ~/ -maxdepth 1 -mindepth 1 -type f -name "*.conf" | sed -e "s/.conf$//")

read -a chipseq_job_list <<< $(qstat -r | grep -Pio "Full\s+jobname:.*_CHIPSEQ$" | sed -e "s/\s//g" | sed -e "s/Fulljobname://" | sed -e "s/_CHIPSEQ$//")

read -a files_to_move <<< $(echo ${chipseq_job_list[@]} | perl -e 'use File::Basename; my $job_info={}; my $conf_info={}; while( my $line = <STDIN> ) { chomp $line; my @array = split( " ", $line ); foreach my $job( @array ) { $job =~ s/^\s+//; $job =~ s/\s+$//; $$job_info{ $job }++; } } foreach my $file( @ARGV ) { chomp $file; $file =~ s/^\s+//; $file =~ s/\s+$//; $$conf_info{ $file }++; } my @files = (); foreach my $file( keys %$conf_info ) { my $base_file = basename( $file ); unless( exists $$job_info{ $base_file } ) { push @files, $file; }} print join( " " , @files );' $(echo ${conf_file_list[@]}))

for file in ${files_to_move[@]}; do
	token=$(basename $file)
	file_path=$(find /zfs/cores/mbcf/mbcf-storage/devel/umv/chipseq/*/*/logs -type f -name "$token.bash") 
	if [ -e "$file_path" ]; then
		copy_to_dir=$(dirname $file_path)
		echo "INFO: Moving chipseq conf file '${file}.conf' to directory '${copy_to_dir}/'" >&2
		cp -f "${file}.conf" "${copy_to_dir}/"
		if [ $? == 0 ]; then
			rm "${file}.conf"
		fi
	fi
done

