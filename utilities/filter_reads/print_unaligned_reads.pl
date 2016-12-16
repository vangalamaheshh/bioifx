#!/usr/bin/perl

# Takes bowtie output and fastq input file
# Outputs fasta file with reads that are not aligned by bowtie
# @author Mahesh Vangala
# @contact vangalamaheshh@gmail.com

use strict;
use warnings;
use Getopt::Long;

my $options = parse_options( );
my $aligned_reads_info = get_aligned_reads_info( $$options{ 'aligned_reads_file' } );
print_unaligned_reads( $aligned_reads_info, $$options{ 'input_reads_file' } );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'aligned_reads_file|a=s', 'input_reads_file|i=s', 'help|h' );
	unless( $$options{ 'aligned_reads_file' } and $$options{ 'input_reads_file' } ) {
		die "Usage: $0 <--aligned_reads_file|-a> <--input_reads_file|-i>\n";
	}
	return $options;
}

sub get_aligned_reads_info {
	my( $file ) = @_;
	my $info = {};
	open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		$line =~ s/\/\d$//;
		$$info{ $line }++;
	}
	close FH or die "Error in clsoign the file, $file, $!\n";
	return $info;
}

sub print_unaligned_reads {
	my( $aligned_reads_info, $input_file ) = @_;
	open( FH, "<$input_file" ) or die "Error in opening file, $input_file, $!\n";
	while( my $header = <FH> ) {
		chomp $header;
		$header =~ s/^@//;
		my $seq = <FH>;
		my $sep = <FH>;
		my $ascii = <FH>;
		unless( exists $$aligned_reads_info{ $header } ) {
			print STDOUT '>' . $header . "\n";
			print STDOUT $seq; 
		}
	}
	close FH or die "Error in closign the file, $input_file, $!\n";
}
