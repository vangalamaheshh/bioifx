#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use XML::Simple;

my $options = parse_options( );
my $run_info_data = XMLin( $$options{ 'run_info_file' }, ForceArray => [ 'Read' ] );
my @indexed_reads = ();
my @nonindexed_reads = ();
foreach my $read( @{ $$run_info_data{ 'Run' }{ 'Reads' }{ 'Read' } } ) {
	if( $$read{ 'IsIndexedRead' } eq 'Y' ) {
		if( $$options{ 'run_type' } eq 'raw_run' ) {
			push @indexed_reads, 'Y'. $$read{ 'NumCycles' };
		}
		else {
			push @indexed_reads, 'I'. $$read{ 'NumCycles' };
		}
	} 
	elsif( $$read{ 'IsIndexedRead' } eq 'N' ) {
		my $num_cycles = int( $$read{ 'NumCycles' } );
		push @nonindexed_reads, 'Y' . ( $num_cycles - 1 ) . 'N';
	}
	else {
		print STDERR "ERROR: There is an anomaly in <Read> in xml file\n";
		exit 1;
	}
}

my $base_mask = '';
if( scalar @nonindexed_reads > 1 ) {
	if( scalar @indexed_reads > 0 ) {
		$base_mask = $nonindexed_reads[ 0 ] . ',' . join( ',', @indexed_reads ) . ',' . $nonindexed_reads[ 1 ];
	}
	else {
		$base_mask = $nonindexed_reads[ 0 ] . ',' . $nonindexed_reads[ 1 ];
	}
}
else {
	if( scalar @indexed_reads > 0 ) {
                $base_mask = $nonindexed_reads[ 0 ] . ',' . join( ',', @indexed_reads );
        }
        else {
                $base_mask = $nonindexed_reads[ 0 ] ;
        }
}

print $base_mask;
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'run_info_file|r=s', 'run_type|t=s', 'help|h' );
	unless( $$options{ 'run_info_file' } and $$options{ 'run_type' } ) {
		die "Usage: $0 <--run_info_file|-r> <--run_type|-t>\n";
	}
	return $options;
}
