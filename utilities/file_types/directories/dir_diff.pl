#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
#--------------------------------------------
# Synopsis: Takes two directory paths.
#           i) Source ii) Destination
#           Outputs the directory names that are present in source path but not in destination path
# 
# @AUTHOR: Mahesh Vangala
#          vangalamaheshh@gmail.com
#          umam_vangala@dfci.harvard.edu
#
# Date: July/07/2014
#--------------------------------------------

my $options = parse_options( );
my $dirs_not_in_dest = get_dirs_not_in_dest( $$options{ 'source_dir' }, $$options{ 'dest_dir' } );
print_diff_dirs( $dirs_not_in_dest );
exit $?;

sub parse_options {
	my $options = {};
	my $usage = "$0 <--source_dir|-s> <--dest_dir|-d>";
	GetOptions( $options, 'source_dir|s=s', 'dest_dir|d=s', 'help|h' );
	unless( $$options{ 'source_dir' } and $$options{ 'dest_dir' } ) {
		die "Usage: $usage\n";
	}
	return $options;
}

sub get_dirs_not_in_dest {
	my( $source_path, $dest_path ) = @_;
	my $diff_dirs = [];
	opendir( my $source_dir, $source_path ) or die "Error in reading the directory, $source_path, $!\n";
	opendir( my $dest_dir, $dest_path ) or die "Error in reading the directory, $dest_path, $!\n";
	my @dest_files = ();
	while( my $dest_file = readdir $dest_dir ) {
		chomp $dest_file;
		push @dest_files, $dest_file;
	} 
	closedir $dest_dir or die "Error in closing the directory, $dest_path, $!\n";
	while( my $cur_dir = readdir( $source_dir ) ) {
		chomp $cur_dir;
		if( ! -e  "$dest_path/$cur_dir" ) { 
			my @files = grep {  /${cur_dir}_.+/ } @dest_files;
			if(  scalar @files == 0 )  {
				push @$diff_dirs, $cur_dir;
			}
		}
	}
	closedir $source_dir or die "Error in closing the directory, $source_path, $!\n";
	return $diff_dirs;
}

sub print_diff_dirs {
	my( $dirs_not_in_dest ) = @_;
	print STDOUT join( " ", @$dirs_not_in_dest );
}
