#!/usr/bin/perl

use strict;
use warnings;

my $matrix_file = shift @ARGV;
my $gtf_file = shift @ARGV;

my ( $header, $matrix_info ) = get_matrix_info( $matrix_file );
my ( $gtf_info ) = get_gtf_info( $gtf_file );
my ( $matrix_with_gene_info ) = get_gene_info( $matrix_info, $gtf_info );
print_matrix_info( $header, $matrix_with_gene_info );

sub get_matrix_info {
	my( $matrix_file ) = @_;
	my $matrix_info = {};
	my $header = undef;
	open( FH, "<$matrix_file" ) or die "Error in opening the file, $matrix_file, $!\n";
	$header = <FH>;
	chomp $header;
	my( $chr_pos, @values ) = split( "\t", $header );
	$header = join( "\t", ( $chr_pos, 'gene_ids', 'transcript_ids', @values ) ); 
	while( my $line = <FH> ) {
		chomp $line;
		my( $chr, @rest ) = split( "\t", $line );
		my( $chr_info, $position ) = split( ":", $chr );
		$$matrix_info{ $chr_info }{ $position }{ 'rest' } = \@rest;
	}
	close FH or die "Error in clsoingt he file, $matrix_file, $!\n";
	return ( $header, $matrix_info );
}

sub get_gtf_info {
	my( $gtf_file ) = @_;
	my $gtf_info = {};
	open( FH, "<$gtf_file" ) or die "Error in opening the file, $gtf_file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my( $chr, $source, $feature, $start, $end, $score, $strand, $frame, $info ) = split( "\t", $line );
		my( $gene_id, $trans_id ) = ( undef, undef );
		if( $info =~ /gene_id\s+\"(.+?)\";.+?transcript_id\s+\"(.+?)\";/ ) {
			( $gene_id, $trans_id ) = ( $1, $2 );
		}
		unless( $gene_id or $trans_id ) {
			die "No Gene_ID and/or Transcript_ID found at line: $.\n";
		}
		push @{ $$gtf_info{ $chr }{ 'start' } }, $start;
		push @{ $$gtf_info{ $chr }{ 'end' } }, $end;
		push @{ $$gtf_info{ $chr }{ 'gene_id' } }, $gene_id;
		push @{ $$gtf_info{ $chr }{ 'trans_id' } }, $trans_id;
	}
	close FH or die "ERror in closign the file, $gtf_file, $!\n";
	return $gtf_info;
}

sub get_gene_info {
	my( $matrix_info, $gtf_info ) = @_;
	my $matrix_with_gene_info = {};
	foreach my $chr( keys %$matrix_info ) {
		print STDERR "INFO: Processing chromosome: $chr\n";
		foreach my $pos( keys %{ $$matrix_info{ $chr } } ) {
			my @indexes = ();
			if( exists $$gtf_info{ $chr } ) {
			INNER_FOR_LOOP:
			foreach my $index( 0 .. scalar @{ $$gtf_info{ $chr }{ 'start' } } - 1 ) {
				if( $pos < ${ $$gtf_info{ $chr }{ 'start' } }[ $index ] ) {
					last INNER_FOR_LOOP;
				}
				if( $pos >= ${ $$gtf_info{ $chr }{ 'start' } }[ $index ] && $pos <= ${ $$gtf_info{ $chr }{ 'end' } }[ $index ] ) {
					push @indexes, $index;
				} 
			}
			}
			my %gene_ids = ();
			my %trans_ids = ();
			unless( @indexes ) {
				%gene_ids = ( '-' => 1 );
				%trans_ids = ( '-' => 1 );
			} else {
				%gene_ids = map { $_ => 1 } @{ $$gtf_info{ $chr }{ 'gene_id' } }[ @indexes ];
				%trans_ids = map { $_ => 1 } @{ $$gtf_info{ $chr }{ 'trans_id' } }[ @indexes ];
			}
			$$matrix_with_gene_info{ $chr . ':' . $pos }{ 'gene_ids' } = join( ',', keys %gene_ids );
			$$matrix_with_gene_info{ $chr . ':' . $pos }{ 'trans_ids' } = join( ',',  keys %trans_ids );
			$$matrix_with_gene_info{ $chr . ':' . $pos }{ 'rest' } = $$matrix_info{ $chr }{ $pos }{ 'rest' };
		}
	}
	return $matrix_with_gene_info;
}

sub print_matrix_info {
	my( $header, $matrix_with_gene_info ) = @_;
	print STDOUT $header, "\n";
	foreach my $pos( sort by_chrom_position keys %$matrix_with_gene_info ) {
		print STDOUT join( "\t", ( $pos, $$matrix_with_gene_info{ $pos }{ 'gene_ids' }, $$matrix_with_gene_info{ $pos }{ 'trans_ids' }, @{ $$matrix_with_gene_info{ $pos }{ 'rest' } } ) ), "\n";
	}	
}

sub by_chrom_position {
        my( $chr_a, $pos_a ) = ( $a =~ /chr(.+)\:(\d+)/ );
        my( $chr_b, $pos_b ) = ( $b =~ /chr(.+)\:(\d+)/ );
        if( $chr_a < $chr_b ) {
                return -1;
        } elsif( $chr_a > $chr_b ) {
                return 1;
        } else {
                return $pos_a <=> $pos_b;
        }
}

