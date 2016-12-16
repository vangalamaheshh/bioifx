#!/usr/bin/perl

use strict;
use warnings;
use lib '/zfs/cores/mbcf/mbcf-storage/devel/umv/software/perl/lib';
use Getopt::Long;
use Text::CSV;

my $options = parse_options( );
generate_raw_sample_sheet( $$options{ 'sample_sheet_in' }, $$options{ 'sample_sheet_out' } );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'sample_sheet_in|i=s', 'sample_sheet_out|o=s', 'help|h' );
	unless( $$options{ 'sample_sheet_in' } and $$options{ 'sample_sheet_out' } ) {
		print "Usage: $0 <--sample_sheet_in|-i> <--sample_sheet_out|-o>\n";
		exit 1;
	}
	return $options;
}

sub generate_raw_sample_sheet {
	my( $in_file, $out_file ) = @_;
	open( my $ifh, "<$in_file" ) or die "Error opening file, $in_file, $!\n";
	open( my $ofh, ">$out_file" ) or die "Error writing to file, $out_file, $!\n";
	
	my $csv = Text::CSV -> new( {
		binary => 1,
		auto_diag => 1,
	} );

	my $header = undef;
	my $index = -1;
	
	WHILE_LOOP:
	while( my $line = $csv -> getline( $ifh ) ) {
		my $token = join( " ", @$line );
		next unless $token;
		if( $token =~ /\bindex\b/i ) {
			$header = $line;
			last WHILE_LOOP;
		}
	}
	foreach my $token( 0 .. scalar @$header - 1 ) {
		if( $$header[ $token ] =~ /index/i ) {
			$index = $token;
		}
	}
	die "Index not found in header" unless $index >= 0;
	print $ofh join( ',', @$header ), "\n";
	my $row = $csv -> getline( $ifh );
	$$row[ 2 ] = $1 . '_' . $2 if( $$row[ 2 ] =~ /(.+?)_.+_(.+)/ );
	print $ofh join( ',', @$row[ 0 .. $index - 1 ] ) . ',,' . join( ',', @$row[ $index + 1 .. scalar @$row - 1 ] ) . "\n";

	close $ofh or die "Error in closign the file, $out_file, $!\n";
	close $ifh or die "Error closing the file, $in_file, $!\n";
}
