#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Spreadsheet::ParseExcel;
use File::Basename;

my $options = parse_options( );
my( $meta_info, $file_info, $source_dir ) = get_meta_info( $$options{ 'metadata_xls_file' } );
my $dest_dir  = "/zfs/cores/mbcf/mbcf-storage/devel/umv/rnaseq";
umask( 0000 );
my( $dir_token ) = ( $$meta_info{ 'run_id' } =~ /^(\d+?)_/ );
$dir_token .= '_' . $$meta_info{ 'ilab_id' };
mkdir( $dest_dir . '/' . $dir_token );
mkdir( $dest_dir . '/' . $dir_token . '/input' );
mkdir( $dest_dir . '/' . $dir_token . '/' . $dir_token . '_rnaseq' );
mkdir( $dest_dir . '/' . $dir_token . '/' . $dir_token . '_rnaseq' . '/logs' );
my $log_dir = $dest_dir . '/' . $dir_token . '/' . $dir_token . '_rnaseq' . '/logs';
my $out_dir = $dest_dir . '/' . $dir_token . '/' . $dir_token . '_rnaseq';
$dest_dir = $dest_dir . '/' . $dir_token . '/input';
copy_input_data( $source_dir, $dest_dir, $file_info );
write_rnaseq_commands( $meta_info, $out_dir, $dest_dir, $log_dir );
execute_bash_scripts( $log_dir );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'metadata_xls_file|m=s', 'help|h' );
	unless( $$options{ 'metadata_xls_file' }  ) {
		print STDERR "Usage: $0 <--metadata_xls_file|-m>\n";
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
	
			elsif( scalar @cols - 1  > 4 && $cols[ 0 ] !~ /sample.*id/i ) {
				my( $sample_id, $sample_name, $left_mate, $sample_type, $grp_num, $grp ) = @cols;
				$grp =~ s/\s+//g;
				$left_mate = `find ${source_dir}/ -name "${left_mate}*_R1*fastq.gz"`;
				chomp $left_mate;
				$left_mate = basename( $left_mate );
				my $right_mate = undef;
				$left_mate =~ s/\.gz$//;
				$$file_info{ $left_mate }++;
				unless( $SE ) {
					$right_mate = $left_mate;
					$right_mate =~ s/_R1/_R2/;
					$$file_info{ $right_mate }++;
				}
				if( $sample_type &&  $sample_type =~ /cont/i ) {
					if( $SE ) {
						push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'single_mates' } }, $left_mate;
					}
					else {
						push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'left_mates' } }, $left_mate;
						push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'right_mates' } }, $right_mate;
					}
					push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'labels' } }, $sample_name;
				}
				elsif( $sample_type =~ /treat/i ) {
					if( $SE ) {
                                        	push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'single_mates' } }, $left_mate;
					}
					else {
						push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'left_mates' } }, $left_mate;
						push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'right_mates' } }, $right_mate;
					}
                                        push @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'labels' } }, $sample_name;
                                }
				else {
					#need to debug this more
					#print STDERR "Only control and treatment types are allowed. $sample_type is not supported.\n";
					#exit( 1 );
				}

			}
			
        	}
	}
	return ( $meta_info, $file_info, $source_dir );
}

sub get_source_dir {
	my( $run_id ) = @_;
	if( $run_id =~ /NS500/ ) {
		return "/zfs/cores/mbcf/mbcf-storage/NS500/analysis/$run_id/concat_per_sample_fastq";
	}
	else {
		if( -e "/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/$run_id/demultiplex_run" ) {
			return "/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/$run_id/demultiplex_run";
		}
		else {
			return "/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/$run_id/concat_per_sample_fastq"
		}
	}
}

sub copy_input_data {
	my( $source_dir, $dest_dir, $file_info ) = @_;
	my $counter = 0;
	my $tracked = {};
	my @fastq_gz_files = split( "\n", `find ${source_dir}/ -name "*.fastq.gz"`);
	foreach my $file( @fastq_gz_files ) {
		my $out_file = basename( $file );
		$out_file =~ s/\.gz$//;
		if( exists $$file_info{ $out_file } and not exists $$tracked{ $out_file } ) {
			$$tracked{ $out_file }++;
			$counter++;
			my $command = "gunzip -c $file 1>$dest_dir/$out_file";
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

sub write_rnaseq_commands {
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
		if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'control' } ) {
			if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'single_mates' } ) {
				@$single_mates_controls = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'single_mates' } };
				foreach my $index( 0 .. scalar @$single_mates_controls -1 ) {
					$control_files[ $index ] = $in_dir . '/' . $$single_mates_controls[ $index ];
				}
			}
			else {
				@$left_mates_controls = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'left_mates' } };
				@$right_mates_controls = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'files' }{ 'right_mates' } };
				foreach my $index( 0 .. scalar @$left_mates_controls - 1 ) {
                        		$control_files[ $index ] = $in_dir . '/' . $$left_mates_controls[ $index ] . '\;'
                                                        . $in_dir . '/' . $$right_mates_controls[ $index ];
                		}	
			}
			@control_labels = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'control' }{ 'labels' } };
			$control_command_string = "--ctrl " . join( ",", @control_files ) . " --ctrllab " . join( ",", @control_labels );
		}
		if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'treat' } ) {
                	if( exists $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'single_mates' } ) {
           	        	@$single_mates_treat = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'single_mates' } };
				foreach my $index( 0 .. scalar @$single_mates_treat - 1 ) {
		                        $treat_files[ $index ] = $in_dir . '/' . $$single_mates_treat[ $index ];
                		}
               		}
                	else {
 	                  	@$left_mates_treat = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'left_mates' } };
        	                @$right_mates_treat = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'files' }{ 'right_mates' } };
				foreach my $index( 0 .. scalar @$left_mates_treat - 1 ) { 
                        		$treat_files[ $index ] = $in_dir . '/' . $$left_mates_treat[ $index ] . '\;'
                                                        . $in_dir . '/' . $$right_mates_treat[ $index ];
                		}
        	        }
			@treat_labels = @{ $$meta_info{ 'grp_name' }{ $grp }{ 'treat' }{ 'labels' } };
			$treat_command_string = "--treat " . join( ",", @treat_files ) . " --treatlab " . join( ",", @treat_labels );
                }

		my $command = "/zfs/cores/mbcf/mbcf-storage/devel/umv/ROOT/bioifx/pipelines/rna_seq/RNA_QC2.py simple " . $control_command_string . " ". $treat_command_string . " -s " 
				. $$meta_info{ 'ref' } . " -n " . $out_dir . '/' . $grp;
		my $out_script = $log_dir . '/' . $grp . '.bash';
		open( OFH, ">$out_script" ) or die "Error writing to the file, $out_script, $!\n";
		print OFH '#!/bin/bash', "\n";
		print OFH '#$ -N ' . $grp . '_RNASEQ' . "\n";
		print OFH '#$ -e ' . $log_dir . '/' . $grp . '.log' . "\n";
		print OFH '#$ -o ' . $log_dir . '/' . $grp . '.out' . "\n";
		print OFH 'echo "Start Time: $(date)" >&2' . "\n"; 
		print OFH 'source /ifs/rcgroups/Apps/RNAseq_pipeline/RNAseqENV/bashrc' . "\n";
		print OFH $command, "\n";
		print OFH 'echo "All analyses done" >&2' . "\n";
		print OFH 'echo "End Time: $(date)" >&2' . "\n";
		close OFH or die "Error closing file, $out_script, $!\n"; 
	}
}

sub execute_bash_scripts {
	my( $script_dir ) = @_;
	opendir(my $dh, $script_dir) or die "can't opendir $script_dir, $!\n";
    	my @bash_scripts = grep { /\.bash$/ && -f "$script_dir/$_" } readdir($dh);
    	closedir $dh;
	foreach my $script( @bash_scripts ) {
		$script = $script_dir . '/' . $script;
		my $command =  "qsub -pe pvm 8 -l mem_free=45G -q all.q $script";
		system( $command );
		unless( $? == 0 ) {
			print STDERR "$command failed with an exit code $?\n";
		}
	}
}
