#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Spreadsheet::ParseExcel;
use File::Basename;

my $options = parse_options( );
my( $meta_info, $file_info, $source_dir ) = get_meta_info( $$options{ 'metadata_xls_file' } );
my $dest_dir  = "/zfs/cores/mbcf/mbcf-storage/devel/umv/chipseq";
umask( 0000 );
my( $dir_token ) = ( $$meta_info{ 'run_id' } =~ /^(\d+?)_/ );
$dir_token .= '_' . $$meta_info{ 'ilab_id' };
mkdir( $dest_dir . '/' . $dir_token );
mkdir( $dest_dir . '/' . $dir_token . '/input' );
mkdir( $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin' );
mkdir( $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin' . '/logs' );
my $log_dir = $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin' . '/logs';
my $out_dir = $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin';
if( $$options{ 'trim' } ) {
	$dest_dir = $out_dir . '/' . 'trimmomatic';
} else {
	$dest_dir = $dest_dir . '/' . $dir_token . '/input';
	copy_input_data( $source_dir, $dest_dir, $file_info );
}
write_chipseq_commands( $meta_info, $out_dir, $dest_dir, $log_dir );
execute_bash_scripts( $log_dir );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'metadata_xls_file|m=s', 'trim|t', 'help|h' );
	unless( $$options{ 'metadata_xls_file' }  ) {
		print STDERR "Usage: $0 <--metadata_xls_file|-m> [--trim|-t]\n";
		exit( 1 );
	}
	if( $$options{ 'metadata_xls_file' } !~ /.*\.xls$/ ) {
		print STDERR "Expected a Micosoft .xls file. But given file is of different format. Exiting ...\n";
		exit( 1 );
	}
	return $options;
}

sub get_meta_info {
	my( $meta_file ) = @_;
	my $meta_info = {};
	my $file_info = {};
	my $SE = undef;
	my $source_dir = undef;
	my $trim_single_token = '.trim.singlemate.fastq.gz';
	my $trim_left_token = '.trim.leftmate.paired.fastq.gz';
	my $trim_right_token = '.trim.rightmate.paired.fastq.gz';
	my $parser = Spreadsheet::ParseExcel -> new( );
	my $workbook = $parser -> parse( $meta_file );
	if( not defined $workbook ) {
        	die $parser -> error( ), "\n";
	}
	foreach my $worksheet( $workbook->worksheets() ) {
		my ( $row_min, $row_max ) = $worksheet->row_range();
        	my ( $col_min, $col_max ) = $worksheet->col_range();

        	for my $row ( $row_min .. $row_max ) {
			my @cols = ();
            		for my $col ( $col_min .. $col_max ) {

                		my $cell = $worksheet->get_cell( $row, $col );
                		next unless $cell;
				
            			push @cols, $cell -> value();
			}
				
			next unless @cols;
	
			if( $cols[ 0 ] =~ /ilab.*id/i ) {
				my $ilab_id = $cols[ 1 ];
				$ilab_id =~ s/^\s+//;
				$ilab_id =~ s/\s+$//;
				$$meta_info{ 'ilab_id' } = $ilab_id;
			}
			
			elsif( $cols[ 0 ] =~ /run.*id/i ) {
				my $run_id = $cols[ 1 ];
				$run_id =~ s/^\s+//;
				$run_id =~ s/\s+$//;
				$$meta_info{ 'run_id' } = $run_id;
				$source_dir = get_source_dir( $$meta_info{ 'run_id' } );	
				my $mate_files = `find ${source_dir}/ -name "*_R2*fastq.gz" | wc -l`;
				$SE = ( $mate_files == 0 ) ? 1 : 0;
			}

			elsif( $cols[ 0 ] =~ /reference/i ) {
				my $ref = $cols[ 1 ];
				$ref =~ s/^\s+//;
				$ref =~ s/\s+$//;
				$$meta_info{ 'ref' } = $ref;
			}

			elsif( $cols[ 0 ] =~ /r.*type/i ) {
				my $run_type = undef;
				my $col_1 = $cols[ 1 ];
				if( $col_1 =~ /tf/i ) {
					$run_type = 'tf';
				} elsif( $col_1 =~ /hist/i ) {
					$run_type = 'histone';
				} elsif( $col_1 =~ /dna/i ) {
					$run_type = 'dnase';
				} else {
					$run_type = 'tf';
				}
				$$meta_info{ 'run_type' } = $run_type;
			}
	
			elsif( scalar @cols - 1  > 4 && $cols[ 0 ] !~ /sample.*id/i ) {
				my( $sample_id, $sample_name, $left_mate, $sample_type, $grp_num, $grp ) = @cols;
				$left_mate =~ s/^\s+//;
				$left_mate =~ s/\s+$//;
				$left_mate = `find ${source_dir}/  -name "${left_mate}*_R1*fastq.gz"`;
				chomp $left_mate;
				$left_mate = basename( $left_mate );
				my $right_mate = undef;
				$$file_info{ $left_mate }++;
				my $trim_file_prefix = $left_mate;
				$trim_file_prefix =~ s/_.*//;
				unless( $SE ) {
					$right_mate = $left_mate;
					$right_mate =~ s/_R1/_R2/;
					$$file_info{ $right_mate }++;
				}
				foreach my $grp_name( split( ',', $grp ) ) {
					if( $sample_type &&  $sample_type =~ /cont/i ) {
						if( $SE ) {
							if( $$options{ 'trim' } ) {
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'control' }{ 'files' }{ 'single_mates' } }, $trim_file_prefix . $trim_single_token;
							} else {
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'control' }{ 'files' }{ 'single_mates' } }, $left_mate;
							}
						}
						else {
							if( $$options{ 'trim' } ) {
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'control' }{ 'files' }{ 'left_mates' } }, $trim_file_prefix . $trim_left_token;
                                                                push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'control' }{ 'files' }{ 'right_mates' } }, $trim_file_prefix . $trim_right_token;

							} else {
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'control' }{ 'files' }{ 'left_mates' } }, $left_mate;
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'control' }{ 'files' }{ 'right_mates' } }, $right_mate;
							}
						}
						push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'control' }{ 'labels' } }, $sample_name;
					}
					elsif( $sample_type =~ /treat/i ) {
						if( $SE ) {
							if( $$options{ 'trim' } ) {
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'treat' }{ 'files' }{ 'single_mates' } }, $trim_file_prefix . $trim_single_token;
							} else {
                                        			push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'treat' }{ 'files' }{ 'single_mates' } }, $left_mate;
							}
						}
						else {
							if( $$options{ 'trim' } ) {
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'treat' }{ 'files' }{ 'left_mates' } }, $trim_file_prefix . $trim_left_token;
                                                                push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'treat' }{ 'files' }{ 'right_mates' } }, $trim_file_prefix . $trim_right_token;

							} else {
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'treat' }{ 'files' }{ 'left_mates' } }, $left_mate;
								push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'treat' }{ 'files' }{ 'right_mates' } }, $right_mate;
							}
						}	
                                        	push @{ $$meta_info{ 'grp_name' }{ $grp_name }{ 'treat' }{ 'labels' } }, $sample_name;
                                	}
					else {
						#need to debug this more
						#print STDERR "Only control and treatment types are allowed. $sample_type is not supported.\n";
						#exit( 1 );
					}
				}

			}
			
        	}
	}
	return ( $meta_info, $file_info, $source_dir );
}

sub get_source_dir {
	my( $run_id ) = @_;
	if( $run_id =~ /(N.50\d)/ ) {
		return "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/$run_id/concat_per_sample_fastq";
	}
	else {
		if( -e "/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/$run_id/demultiplex_run" ) {
			return "/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/$run_id/demultiplex_run";
		} else {
			return "/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/$run_id/concat_per_sample_fastq"
		}
	}
}

sub copy_input_data {
	my( $source_dir, $dest_dir, $file_info ) = @_;
	my $counter = 0;
	my @fastq_gz_files = split( "\n", `find ${source_dir}/ -name "*.fastq.gz"`);
	foreach my $file( @fastq_gz_files ) {
		my $out_file = basename( $file );
		#$out_file =~ s/\.gz$//;
		if( exists $$file_info{ $out_file } ) {
			$counter++;
			#my $command = "gunzip -c $file 1>$dest_dir/$out_file";
			my $command = "ln -s $file ${dest_dir}/${out_file}";
			system( $command );
			unless( $? == 0 ) {
				print STDERR "$command failed with exit code: $?\n";
				exit( 1 );
			}
		}
	}	
	
	if( $counter != scalar keys %$file_info ) {
		print STDERR "All files mentioned in the metadata sheet won't match with files found. Exiting ...\n";
		exit( 1 );
	}
}

sub write_chipseq_commands {
	my( $meta_info, $out_dir, $in_dir, $log_dir ) = @_;
        foreach my $grp( keys %{ $$meta_info{ 'grp_name' } } ) {
                my( $single_mates_controls, $left_mates_controls, $right_mates_controls ) = ( [], [], [] );
                my( $single_mates_treat, $left_mates_treat, $right_mates_treat ) = ( [], [], [] );
                my @control_files = ();
                my @control_labels = ();
                my @treat_files = ();
                my @treat_labels = ();
                my $control_command_string = "";
                my $treat_command_string = "";
		my $paired_end = 0;
                if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'control' } ) {
                        if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'single_mates' } ) {
                                @$single_mates_controls = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'single_mates' } };
                                foreach my $index( 0 .. scalar @$single_mates_controls -1 ) {
                                        $control_files[ $index ] = $in_dir . '/' . $$single_mates_controls[ $index ];
                                }
                        }
                        else {  
				$paired_end = 1;
                                @$left_mates_controls = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'left_mates' } };
                                @$right_mates_controls = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'right_mates' } };
                                foreach my $index( 0 .. scalar @$left_mates_controls - 1 ) {
                                        $control_files[ $index ] = $in_dir . '/' . $$left_mates_controls[ $index ] . ','
                                                        . $in_dir . '/' . $$right_mates_controls[ $index ];
                                }
                        }
                        @control_labels = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'labels' } };
                        $control_command_string = "--cont " . join( "\;", @control_files );
                }
                if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'treat' } ) {
                        if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'single_mates' } ) {
                                @$single_mates_treat = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'single_mates' } };
                                foreach my $index( 0 .. scalar @$single_mates_treat - 1 ) {
                                        $treat_files[ $index ] = $in_dir . '/' . $$single_mates_treat[ $index ];
                                }
                        }
                        else {  
				$paired_end = 1;
                                @$left_mates_treat = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'left_mates' } };
                                @$right_mates_treat = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'right_mates' } };
                                foreach my $index( 0 .. scalar @$left_mates_treat - 1 ) {
                                        $treat_files[ $index ] = $in_dir . '/' . $$left_mates_treat[ $index ] . ','
                                                        . $in_dir . '/' . $$right_mates_treat[ $index ];
                                }
                        }
                        @treat_labels = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'labels' } };
                        $treat_command_string = "--treat " . join( "\;", @treat_files );
                }
		my $paired_end_token = $paired_end ? '--pe' : '';
		unless( exists $$meta_info{ 'run_type' } ) {
			$$meta_info{ 'run_type' } = 'tf';
		}
		my $command = "ChiLin2 simple $paired_end_token -p both --threads 8  $control_command_string $treat_command_string --rtype " . 
			$$meta_info{ 'run_type' }   . " -i " . $grp . " -s " . $$meta_info{ 'ref' } . " -o " . $out_dir . '/' . $grp;
		my $out_script = $log_dir . '/' . $grp . '.chipseq.bash';
		open( OFH, ">$out_script" ) or die "Error writing to the file, $out_script, $!\n";
		print OFH '#!/bin/bash', "\n";
		print OFH '#$ -N ' . $grp . '_CHIPSEQ' . "\n";
		print OFH '#$ -e ' . $log_dir . '/' . $grp . '.log' . "\n";
		print OFH '#$ -o ' . $log_dir . '/' . $grp . '.out' . "\n";
		print OFH 'echo "Start Time: $(date)" >&2' . "\n";
		print OFH 'hostname >&2' . "\n"; 
		print OFH 'source /ifs/rcgroups/Apps/env_chilin2/bin/activate' . "\n";
		print OFH 'export PATH=/apps/jdk1.8.0_20/bin:/apps/R-2.14.1c/bin:$PATH' . "\n";
		print OFH $command, "\n";
		print OFH 'echo "All analyses done" >&2' . "\n";
		print OFH 'echo "End Time: $(date)" >&2' . "\n";
		close OFH or die "Error closing file, $out_script, $!\n"; 
	}
}

sub execute_bash_scripts {
	my( $script_dir ) = @_;
	opendir(my $dh, $script_dir) or die "can't opendir $script_dir, $!\n";
    	my @bash_scripts = grep { /\.chipseq\.bash$/ && -f "$script_dir/$_" } readdir($dh);
    	closedir $dh;
	foreach my $script( @bash_scripts ) {
		$script = $script_dir . '/' . $script;
		my $command =  "qsub -pe pvm 8 -l mem_free=35G -q all.q $script";
		system( $command );
		unless( $? == 0 ) {
			print STDERR "$command failed with an exit code $?\n";
		}
	}
}
