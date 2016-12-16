#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

print_info( get_info( ) );
exit $?;

sub get_info {
	my $info = {};
	my @samples = ();
	while( my $file = <STDIN> ) {
		chomp $file;
		open( VCF, "<$file" ) or die "Error in opening the file, $file, $!\n";
		my $sample_name = basename( dirname( $file ) );
		while( my $line = <VCF> ) {
			chomp $line;
			next if $line =~ /^#/;
			my( $chrom, $pos, $id, $ref, $alt, $qual ) = split( "\t", $line );
			$alt =~ s/,/;/g;
			$id =~ s/,/;/g;
			$$info{ $chrom . ':' . $pos }{ 'ref' } = $ref;
			$$info{ $chrom . ':' . $pos }{ 'snp_id' } = $id;
			$$info{ $chrom . ':' . $pos }{ 'sample_name' }{ $sample_name } = { 'alt' => $alt, 'qual' => $qual };
		}
		
		push @samples, $sample_name;
		close VCF or die "Error in clsoing the file, $file, $!\n";
	}
	return $info, @samples;
}

sub print_info {
	my( $info, @samples ) = @_;
	print STDOUT join( ",", ( "Position", "Reference", "SNP_ID", @samples ) ), "\n";
	foreach my $pos( sort by_chrom_position keys %$info ) {
		my @values = ();
		foreach my $sample( @samples ) {
			unless( exists $$info{ $pos }{ 'sample_name' }{ $sample } ) {
				push @values, '-';
			} else {
				push @values, $$info{ $pos }{ 'sample_name' }{ $sample }{ 'alt' } . ' (' . $$info{ $pos }{ 'sample_name' }{ $sample }{ 'qual' } . ')';  
			}
		}
		print STDOUT join( ",", ( $pos, $$info{ $pos }{ 'ref' }, $$info{ $pos }{ 'snp_id' }, @values ) ), "\n";
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
