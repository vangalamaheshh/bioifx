#!/usr/bin/perl

#----------------------------------------------------
# Description:
# ------------
# Illumina nextseq generates 4 fastq files, one file per lane, for each sample.
# This script outputs one fastq file per sample by concatenating all lane files
# for that sample.
# @AUTHOR: Mahesh Vangala
# @email : vangalamaheshh@gmail.com
# Date: July 01, 2014
#----------------------------------------------------- 
use strict;
use warnings;
use Getopt::Long;
use Log::Log4perl qw(:easy);
use threads;
use threads::shared;

#------------ Initializing logger ---------------------
Log::Log4perl -> easy_init( $INFO );
#------------------------------------------------------

my $options = parse_options( );
my @fq_files = get_fastq_files( $$options{ 'input_dir' } );
INFO "MSG: Number of fastq.gz files present: ". scalar @fq_files;
my $file_info = get_file_info( @fq_files );
$file_info = concat_files( $file_info, $$options{ 'input_dir' }, $$options{ 'concat_dir' } );
INFO "MSG: Concatenating fastq files is done";
INFO "MSG: Success";
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'input_dir|i=s', 'concat_dir|c=s', 'help|h' );
	unless( $$options{ 'input_dir' } and $$options{ 'concat_dir' } ) {
		die "Usage: $0 <--input_dir|-i> <--concat_dir|-c>\n";
	}
	return $options;
}

sub get_fastq_files {
	my( $input_dir ) = @_;
	opendir( my $dh, $input_dir ) or die "Error in reading directory, $input_dir, $!\n";
	my @fq_files = grep { /.*\.fastq\.gz$/ && -f "$input_dir/$_" } readdir( $dh );
	closedir( $dh ) or die "Error in closing directory, $input_dir, $!\n";
	return @fq_files;
}

sub get_file_info {
	my( @fq_files ) = @_;
	my $file_info = {};
	foreach my $fq_file( @fq_files ) {
		if( $fq_file =~ /(.+)_L00(\d)_(R\d)/ ) {
			my( $sample, $lane, $mate ) = ( $1, $2, $3 );
			$$file_info{ $sample }{ $mate }{ $lane } = $fq_file;	
		}
		elsif( $fq_file =~ /(.+)_L00(\d)_(I1)/ ) {
			my( $sample, $lane, $index ) = ( $1, $2, $3 );
			$$file_info{ $sample }{ $index }{ $lane } = $fq_file;
		}
		elsif( $fq_file =~ /(.+)_L00(\d)_(I2)/ ) {
                        my( $sample, $lane, $index ) = ( $1, $2, $3 );
                        $$file_info{ $sample }{ $index }{ $lane } = $fq_file;
                }
		else {
			INFO "ERROR: $fq_file doesn't follow convention 'samplename_Lane_mate'";
			exit 1;
		}
	}
	return $file_info;
}

sub concat_files {
	my( $file_info, $input_dir, $concat_dir ) = @_;
	my @threads = ();
	foreach my $sample( keys %$file_info ) {
		INFO "MSG: Processing sample: $sample";
		INFO "--------------------------------------";
		INFO "MSG: concatenating 'LeftMate' files";
		my @leftmates = ();
		foreach my $lane( sort { $a <=> $b } keys %{ $$file_info{ $sample }{ 'R1' } } ) {
			push @leftmates, $input_dir. "/" . $$file_info{ $sample }{ 'R1' }{ $lane };
		}
		my $output_file_leftmate = $concat_dir . "/" . $sample . "_" . "R1" . ".fastq.gz";
		my $leftmate_list = join( " ", @leftmates );
		my $command = "zcat $leftmate_list | gzip 1>$output_file_leftmate";
		my $thread = threads -> new( \&execute_command, $command );
		push @threads, $thread;
		
		if( exists $$file_info{ $sample }{ 'R2' } ) {
			my @rightmates = ();
			INFO "MSG: concatenating 'RightMate' files";
			foreach my $lane( sort { $a <=> $b } keys %{ $$file_info{ $sample }{ 'R2' } } ) {
				push @rightmates, $input_dir . "/" . $$file_info{ $sample }{ 'R2' }{ $lane };
			}
			my $output_file_rightmate = $concat_dir . "/" . $sample . "_" . "R2" . ".fastq.gz";
			my $rightmate_list = join( " ", @rightmates );
			my $command = "zcat $rightmate_list | gzip 1>$output_file_rightmate";
			my $thread = threads -> new( \&execute_command, $command );
			push @threads, $thread;
		}

		if( exists $$file_info{ $sample }{ 'I1' } ) {
                        my @index_files = ();
                        INFO "MSG: concatenating 'Index I1' files";
                        foreach my $lane( sort { $a <=> $b } keys %{ $$file_info{ $sample }{ 'I1' } } ) {
                                push @index_files, $input_dir . "/" . $$file_info{ $sample }{ 'I1' }{ $lane };
                        }
                        my $output_file_index = $concat_dir . "/" . $sample . "_" . "I1" . ".fastq.gz";
                        my $index_list = join( " ", @index_files );
                        my $command = "zcat $index_list | gzip 1>$output_file_index";
                        my $thread = threads -> new( \&execute_command, $command );
                        push @threads, $thread;
                }

		if( exists $$file_info{ $sample }{ 'I2' } ) {
                        my @index_files = ();
                        INFO "MSG: concatenating 'Index I2' files";
                        foreach my $lane( sort { $a <=> $b } keys %{ $$file_info{ $sample }{ 'I2' } } ) {
                                push @index_files, $input_dir . "/" . $$file_info{ $sample }{ 'I2' }{ $lane };
                        }
                        my $output_file_index = $concat_dir . "/" . $sample . "_" . "I2" . ".fastq.gz";
                        my $index_list = join( " ", @index_files );
                        my $command = "zcat $index_list | gzip 1>$output_file_index";
                        my $thread = threads -> new( \&execute_command, $command );
                        push @threads, $thread;
                }

	}
	$_ -> join foreach @threads;
	return $file_info;
}

sub execute_command {
	my( $command ) = @_;
	my $exit_code = system( $command );
	if( $exit_code ) {
		INFO "ERROR: $command failed with exit code: $exit_code";
	}
	else {
		INFO "MSG: $command succeeded with exit code '0'";
	}
}
