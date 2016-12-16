#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename; 

my $options = parse_options( );
my ($info,$track) = process_info( $$options{ 'sorted_bam' } );
print_info( $info, $track );
exit $?;

sub parse_options {
	my $options = {};
	GetOptions( $options, 'sorted_bam|b=s@', 'help|h' );
	unless( $$options{ 'sorted_bam' } ) {
		print STDERR "Usage: $0 <--sorted_bam|-b> [--sorted_bam|-b]\n";
		exit 1;
	}
	return $options;
}

sub process_info {
	my( $files ) = @_;
	my $info = {};
	my $track = {};
	foreach my $file( @$files ) {
		open( FH, "samtools view -S $file |" ) or die "Error in opening the file, $file, $!\n";
		my $base_file = basename( $file );
		$base_file =~ s/\.sorted\.bam//;
		while( my $line = <FH> ) {
			chomp $line;
			next if $line =~ /^\@/;
			my( $query_id, $flag, $sub_id ) = split( "\t", $line );
			if( $flag == 16 || $flag == 0 ) {
				$$info{ $base_file }{ $sub_id }++;
				$$track{$sub_id}++;
			}
		}
		close FH or die "Error in closing the file, $file, $!\n";
	}
	return ($info, $track);
}

sub print_info {
	my( $info, $track ) = @_;
	my @files = keys %$info;
	my @base_names = ();
	foreach my $file( @files ) {
		my $base = basename( $file );
		push @base_names, $base;
	}
	print STDOUT join( ",", ( 'miRNA_ID', @base_names ) ), "\n";
	foreach my $miRNA( keys %$track ){
		my @nums = ();
		foreach my $file( @files ) {
			if( exists $$info{ $file }{ $miRNA } ) {
				push @nums, $$info{ $file }{ $miRNA };
			} else {
				push @nums, 0;
			}
		}
		print STDOUT join( ',', ( $miRNA, @nums ) ), "\n";
	}
}
