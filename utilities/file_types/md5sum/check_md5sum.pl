#!/usr/bin/perl

#---------------------------------
# Description: Takes two files with md5sum info.
#              Checks for equality of md5sum info.
#              Outputs the filenames that have diff md5sum
#              using token MD5SUM_DIFFERENCE:
#              and the files that are not found using token
#              FILE_NOT_FOUND:
# @AUTHOR: Mahesh Vangala
#          vangalamaheshh@gmail.com
#          umam_vangala@dfci.harvard.edu
# @DATE: July, 09, 2014
#-----------------------------------

use strict;
use warnings;
use Getopt::Long;

my $options = parse_options( );
my $mbcftp_info = get_mbcftp_info( $$options{ 'mbcftp_checksum_file' } );
my $rchpc1_info = get_rchpc1_info( $$options{ 'rchpc1_checksum_file' } );
print_diffs( $mbcftp_info, $rchpc1_info );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'mbcftp_checksum_file|m=s', 'rchpc1_checksum_file|r=s', 'help|h' );
	my $usage = "$0 <--mbcftp_checksum_file|-m> <--rchpc1_checksum_file|-r>";
	unless( $$options{ 'mbcftp_checksum_file' } and $$options{ 'rchpc1_checksum_file' } ) {
		die "Usage: $usage\n";
	}
	return $options;
}

sub get_mbcftp_info {
	my( $file ) = @_;
	my $info = {};
	open( my $fh, "<$file" ) or die "Error in opening the file, $file, $!\n";
	while( my $line = <$fh> ) {
		chomp $line;
		if( $line =~ /MD5\s+\((.+)\)\s+=\s+(\S+)/ ) {
			my( $file_name, $md5sum ) = ( $1, $2 );
			next if $file_name =~ /\.DS_Store$/;
			$$info{ $file_name }{ 'tracked' } = 0;
			$$info{ $file_name }{ 'md5sum' } = $md5sum;
		}
	}
	close $fh or die "Error in clsoign the file, $file, $!\n";
	return $info;
}

sub get_rchpc1_info {
	my( $file ) = @_;
	my $info = {};
	open( my $fh, "<$file" ) or die "Error in opening the file, $file, $!\n";
	while( my $line = <$fh> ) {
		chomp $line;
		if( $line =~ /(\S+)\s+(.+)/ ) {
			my( $file_name, $md5sum ) = ( $2, $1 );
			next if $file_name =~ /\.DS_Store$/;
			$$info{ $file_name }{ 'tracked' } = 0;
			$$info{ $file_name}{ 'md5sum' } = $md5sum;
		}
	}
	close $fh or die "Error in clsing the file, $file, $!\n";
	return $info;
}


sub print_diffs {
	my( $mbcftp_info, $rchpc1_info ) = @_;
	foreach my $file_name( keys %$mbcftp_info ) {
		if( exists $$rchpc1_info{ $file_name } ) {
			$$mbcftp_info{ $file_name }{ 'tracked' } = 1;
			unless( $$mbcftp_info{ $file_name }{ 'md5sum' } eq
				$$rchpc1_info{ $file_name }{ 'md5sum' } ) {
				print STDOUT "MD5_DIFFERENCE: $file_name\n";
			}
		}
		else {
			print STDOUT "FILE_NOT_FOUND: $file_name\n";
		}
	}
}
