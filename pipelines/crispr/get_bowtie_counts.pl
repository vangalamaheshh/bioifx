#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename; 

my $options = parse_options( );
my $info = process_info( $$options{ 'bowtie_out_file' } );
print_info( $info, $$options{ 'ref_input_file' } );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'bowtie_out_file|b=s@', 'ref_input_file|r=s', 'help|h' );
	unless( $$options{ 'bowtie_out_file' } and $$options{ 'ref_input_file' } ) {
		print STDERR "Usage: $0 <--bowtie_out_file|-b> <--ref_input_file|-r> [--bowtie_out_file|-b]\n";
		exit 1;
	}
	return $options;
}

sub process_info {
	my( $files ) = @_;
	my $info = {};
	foreach my $file( @$files ) {
		open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
		while( my $line = <FH> ) {
			chomp $line;
			next if $line =~ /^\@/;
			my( $query_id, $flag, $sub_id ) = split( "\t", $line );
			next if $sub_id eq '*';
			my $base_file = basename( $file );
			$$info{ $base_file }{ $sub_id }++;
		}
		close FH or die "Error in closing the file, $file, $!\n";
	}
	return $info;
}

sub print_info {
	my( $info, $ref_file ) = @_;
	my @files = keys %$info;
	my @base_names = ();
	foreach my $file( @files ) {
		push @base_names, basename( $file );
	}
	print STDOUT join( ",", ( 'Gene', 'Symbol', 'OligoSeq', @base_names ) ), "\n";
	open( FH, "<$ref_file" ) or die "Error in opening the file, $ref_file, $!\n";
	while( my $line = <FH> ) {
		next unless substr( $line, 0, 1 ) eq '>' ;
		chomp $line;
		$line =~ s/^>//;
		my @nums = ();
		foreach my $file( @files ) {
			if( exists $$info{ $file }{ $line } ) {
				push @nums, $$info{ $file }{ $line };
			} else {
				push @nums, 0;
			}
		}
		print STDOUT join( ',', ( split( ";", $line ), @nums ) ), "\n";
	}
	close FH or die "Error in closing the file, $ref_file, $!\n";
}
