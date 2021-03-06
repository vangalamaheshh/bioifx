#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Spreadsheet::ParseExcel;
use File::Basename;

my $options = parse_options( );
my $trim_info = get_trim_info( $$options{ 'metadata_xls_file' } );
my $dest_dir  = "/zfs/cores/mbcf/mbcf-storage/devel/umv/chipseq";
umask( 0000 );
my( $dir_token ) = ( $$trim_info{ 'run_id' } =~ /^(\d+?)_/ );
$dir_token .= '_' . $$trim_info{ 'ilab_id' };
mkdir( $dest_dir . '/' . $dir_token );
mkdir( $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin' );
mkdir( $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin' . '/logs' );
my $log_dir = $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin' . '/logs';
my $out_dir = $dest_dir . '/' . $dir_token . '/' . $dir_token . '_chilin';
my $trim_dir = $out_dir . '/' . 'trimmomatic';
mkdir( $trim_dir );
write_trimmomatic_commands( $trim_info, $trim_dir, $log_dir );
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

sub get_trim_info {
	my( $meta_file ) = @_;
	my $trim_info = {};
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
				$$trim_info{ 'ilab_id' } = $ilab_id;
			}
			
			elsif( $cols[ 0 ] =~ /run.*id/i ) {
				my $run_id = $cols[ 1 ];
				$run_id =~ s/^\s+//;
				$run_id =~ s/\s+$//;
				$$trim_info{ 'run_id' } = $run_id;
				$source_dir = get_source_dir( $$trim_info{ 'run_id' } );	
				my $mate_files = `find ${source_dir}/ -name "*_R2*fastq.gz" | wc -l`;
				$SE = ( $mate_files == 0 ) ? 1 : 0;
			}
	
			elsif( scalar @cols - 1  > 4 && $cols[ 0 ] !~ /sample.*id/i ) {
				my( $sample_id, $sample_name, $left_mate, $sample_type, $grp_num, $grp ) = @cols;
				$left_mate =~ s/^\s+//;
				$left_mate =~ s/\s+$//;
				$left_mate = `find ${source_dir}/  -name "${left_mate}*_R1*fastq.gz"`;
				chomp $left_mate;
				my $dir_name = dirname( $left_mate );
				$left_mate = basename( $left_mate );
				my $right_mate = undef;
				my $trim_file_prefix = $left_mate;
				$trim_file_prefix =~ s/_.*//;
				$$trim_info{ 'trim_in_dir' } = $dir_name;
				$$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'count' }++;	
				$$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'left_mate' } = $left_mate;
				unless( $SE ) {
					$right_mate = $left_mate;
					$right_mate =~ s/_R1/_R2/;
					$$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'PE' } = 1;
					$$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'right_mate' } = $right_mate;
				}
			} else {
#				print STDERR join( "\t", @cols ), "\n";
			}
			
        	}
	}
	return $trim_info;
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
			return "/zfs/cores/mbcf/mbcf-storage/MiSeq/analysis/$run_id/concat_per_sample_fastq";
		}
	}
}

sub write_trimmomatic_commands {
	my( $trim_info, $trim_out_dir, $log_dir ) = @_;
	my $java_exec = '/ifs/rcgroups/mbcf/umv/software/jdk1.8.0_45/bin/java';
	my $trim_exec = '/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/trimmomatic-0.32.jar';
      	my $pe_adapter_file='/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TRUSEQ_NEXTERA_PE.fa';
	my $se_adapter_file = '/zfs/cores/mbcf/mbcf-storage/devel/umv/software/Trimmomatic-0.32/adapters/TRUSEQ_NEXTERA_SE.fa';
	my $trim_in_dir = $$trim_info{ 'trim_in_dir' };

	my $job_count = 1;
        foreach my $trim_file_prefix( keys %{ $$trim_info{ 'trim_files' } } ) {
		my $command = undef;
		if( exists $$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'PE' } ) {
			my $left_mate_in = $trim_in_dir . '/' . $$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'left_mate' };
			my $right_mate_in = $trim_in_dir . '/' . $$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'right_mate' };
			my $left_mate_paired_out = $trim_out_dir . '/' . $trim_file_prefix . '.trim.leftmate.paired.fastq.gz';
			my $left_mate_unpaired_out = $trim_out_dir . '/' . $trim_file_prefix . '.trim.leftmate.unpaired.fastq.gz';
			my $right_mate_paired_out = $trim_out_dir . '/' . $trim_file_prefix . '.trim.rightmate.paired.fastq.gz';
			my $right_mate_unpaired_out = $trim_out_dir . '/' . $trim_file_prefix . '.trim.rightmate.unpaired.fastq.gz';
			
			$command = "$java_exec -jar $trim_exec PE -phred33 -threads 8 $left_mate_in $right_mate_in $left_mate_paired_out $left_mate_unpaired_out $right_mate_paired_out $right_mate_unpaired_out ILLUMINACLIP:${pe_adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36";
		}	
		else {
			my $single_mate_in = $trim_in_dir . '/' . $$trim_info{ 'trim_files' }{ $trim_file_prefix }{ 'left_mate' };
			my $single_mate_out = $trim_out_dir . '/' . $trim_file_prefix . '.trim.singlemate.fastq.gz';
			
			$command = "$java_exec -jar $trim_exec SE -phred33 -threads 8 $single_mate_in $single_mate_out ILLUMINACLIP:${se_adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36";
		}
			
		my $out_script = $log_dir . '/' . $trim_file_prefix . '.trim.bash';
		open( OFH, ">$out_script" ) or die "Error writing to the file, $out_script, $!\n";
		print OFH '#!/bin/bash', "\n";
		print OFH '#$ -N ' . 'TRIM_' . $job_count . "\n";
		print OFH '#$ -e ' . $log_dir . '/' . $trim_file_prefix . '.log' . "\n";
		print OFH '#$ -o ' . $log_dir . '/' . $trim_file_prefix . '.out' . "\n";
		print OFH 'echo "Start Time: $(date)" >&2' . "\n";
		print OFH 'hostname >&2' . "\n"; 
		print OFH $command, "\n";
		print OFH 'echo "All analyses done" >&2' . "\n";
		print OFH 'echo "End Time: $(date)" >&2' . "\n";
		close OFH or die "Error closing file, $out_script, $!\n";
		$job_count++; 
	}
}

sub execute_bash_scripts {
	my( $script_dir ) = @_;
	opendir(my $dh, $script_dir) or die "can't opendir $script_dir, $!\n";
    	my @bash_scripts = grep { /\.trim\.bash$/ && -f "$script_dir/$_" } readdir($dh);
    	closedir $dh;
	foreach my $script( @bash_scripts ) {
		$script = $script_dir . '/' . $script;
		my $command =  "qsub -pe pvm 8 -q all.q $script";
		system( $command );
		unless( $? == 0 ) {
			print STDERR "$command failed with an exit code $?\n";
		}
	}
}
