#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

=head1 SYNOPSIS

   Takes the output fasta file from amplicon noise pipeline and 
   generates one fasta file per sample.

   PARAMETERS:
  
   --fasta_file|-f   ==> output fasta file from amplicon noise pipeline
   --out_dir|-o      ==> output directory where you want to write the output files
   --suffix|-s       ==> Optional. Append the fasta header and output filename with this suffix.
   --help|-h         ==> print this help message and exit

=head1 AUTHOR
   
   Mahesh Vangala
   vangalamaheshh@gmail.com

=cut

my $options = parse_options( );
create_dir( $$options{ 'out_dir' } );
my $info = get_info( $$options{ 'fasta_file' }, $$options{ 'suffix' } );
output_files( $info, $$options{ 'out_dir' } );
exit $?;

sub create_dir {
	my( $dir_name ) = @_;
	if( ! -e "$dir_name" ) {
		mkdir "$dir_name";
	}
	elsif( ! -d "$dir_name" ) {
		print STDERR "A file exists with name, $dir_name, $!\n";
	}
}

sub parse_options {
	my $options = {};
	GetOptions( $options, 'fasta_file|f=s', 'out_dir|o=s', 'suffix|s=s', 'help|h' ) ||
		pod2usage( { -output => *STDOUT, -verbose => 2 } );
	unless( $$options{ 'fasta_file' } and $$options{ 'out_dir' } ) {
		pod2usage( );
	}
	
	return $options;
}

sub get_info {
	my( $fasta_file, $suffix ) = @_;
	my $info = {};
	open( FASTA, "<$fasta_file" ) or die "Error opening the file, $fasta_file, $!\n";
	OUTER_WHILE:
	while( my $line = <FASTA> ) {
		my $token = undef;
		if( $line =~ /^>(.+)_/ ) {
			if( $suffix ) {
				$token = "$1_$suffix";
				chomp $line;
				$line .= "_$suffix\n";
			}
			else {
				$token = $1;
			}
		}
		push @{ $$info{ $token } }, $line;
		while( $line = <FASTA> ) {
			if( $line =~ /^>/ ) { redo OUTER_WHILE; }
			push @{ $$info{ $token } }, $line;
		}
		last OUTER_WHILE;
	}
	close FASTA or die "Error opening the file, $fasta_file, $!\n";
	return $info;
}

sub output_files {
	my( $info, $out_dir ) = @_;
	while( my ( $token, $ref_array ) = each %$info ) {
		my $output_file = "$out_dir/$token.fasta";
		open( OFH, ">$output_file" ) or die "Error in writing to the file, $output_file, $!\n";
		print OFH @$ref_array;
		close OFH or die "Error in closing the file, $output_file, $!\n";
	}
}
