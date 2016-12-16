#!/usr/bin/perl
#
use strict;
use warnings;

my $barcode_1 = shift @ARGV;

my $uni_adapter = "TCTTGTGGAAAGGACGAAACACCG";

while( my $header = <STDIN> ) {
	my $seq = <STDIN>;
	my $sep = <STDIN>;
	my $ascii = <STDIN>;
	
	if( $seq =~ /(.+${barcode_1}${uni_adapter})(....................)/ ) {
		print STDOUT $header;
		print STDERR $header;
		print STDOUT $2, "\n";
		print STDERR $seq;
		print STDOUT $sep;
		print STDERR $sep;
		my $start = length( $1 );
		print STDOUT substr( $ascii, $start, 20 ), "\n";
		print STDERR $ascii;
	}
	#else {
	#	print STDERR $header;
	#	print STDERR $seq;
	#	print STDERR $sep;
	#	print STDERR $ascii;
	#}
}


