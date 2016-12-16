#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $options = parse_options();
my $info = get_info_after_removing_dups_and_revComps( $$options{ 'ref_fasta' } );
print_info( $info );
exit $?;

sub print_info {
	my( $info ) = @_;
	foreach my $seq( keys %$info ) {
		print STDOUT join( ';', ( $$info{ $seq }{ 'symbol' }, $$info{ $seq }{ 'id' }, $seq ) ) . "\n";
		print STDOUT $seq, "\n";
	}
}

sub parse_options {
	my $options = {};
	GetOptions( $options, 'ref_fasta|r=s', 'help|h' );
	unless( $$options{ 'ref_fasta' } ) {
		my $usage = "$0 <--ref_fasta|-r>";
		print STDERR $usage, "\n";
		exit 1;
	}
	return $options;
}

sub get_info_after_removing_dups_and_revComps {
	my( $ref_fasta ) = @_;
	my $info = {};
	open( FH, "<$ref_fasta" ) or die "Error in opening the file, $ref_fasta, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my( $symbol, $id ) = split( ";", $line );
		my $seq = <FH>;
		chomp $seq;
		my $rev_comp = $seq;
		$rev_comp =~ tr/AGCT/TCGA/;
		$rev_comp = reverse( $rev_comp );
		if( exists $$info{ $seq } ) {
			remove_dup_info( $info, $symbol, $id, $seq, undef );
		}
		elsif( exists $$info{ $rev_comp } ) {
			remove_dup_info( $info, $symbol, $id, $rev_comp, $seq );
		} 
		else {
			$$info{ $seq }{ 'symbol' } = $symbol;
			$$info{ $seq }{ 'id' } = $id;
		}
	}
	close FH or die "Error in closing the file, $ref_fasta, $!\n";
	return $info;
}

sub remove_dup_info {
	my( $info, $cur_symbol, $cur_id, $ref_seq, $cur_seq ) = @_;
	my( $cur_string, $cur_num ) = ( $cur_id =~ /(\D+)(\d+)/ );
	my( $ref_string, $ref_num ) = ( $$info{ $ref_seq }{ 'id' } =~ /(\D+)(\d+)/ );
	if( $cur_string lt $ref_string ) {
		if( $cur_seq ) {
			print STDERR "REV_COMP: " . join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) .
                " replaced with " . join( ';', ( $cur_symbol, $cur_id, $cur_seq ) ) . "\n";
			delete $$info{ $ref_seq };
			$$info{ $cur_seq }{ 'symbol' } = $cur_symbol;
			$$info{ $cur_seq }{ 'id' } = $cur_id;
		} else {
			print STDERR "DUPS: " . join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) .
				" replaced with " . join( ';', ( $cur_symbol, $cur_id, $ref_seq ) ) . "\n";
			$$info{ $ref_seq }{ 'id' } = $cur_id;
			$$info{ $ref_seq }{ 'symbol' } = $cur_symbol;
		}
	}
	elsif( $ref_string lt $cur_string ) {
		if( $cur_seq ) { 
			print STDERR "REV_COMP: " . join( ';', ( $cur_symbol, $cur_id, $cur_seq ) ) . " replaced with " . 
				join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) . "\n";
		} else {
			print STDERR "DUPS: ". join( ';', ( $cur_symbol, $cur_id, $ref_seq ) ) . 
				" replaced with " . join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) . "\n";
		}
	}
	else {
		if( $cur_num lt $ref_num ) {
			if( $cur_seq ) {
				print STDERR "REV_COMP: " . join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) .
                	" replaced with " . join( ';', ( $cur_symbol, $cur_id, $cur_seq ) ) . "\n";
            	delete $$info{ $ref_seq };
            	$$info{ $cur_seq }{ 'symbol' } = $cur_symbol;
            	$$info{ $cur_seq }{ 'id' } = $cur_id;
			} else {
     	   		print STDERR "DUPS: " . join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) .
        	    	" replaced with " . join( ';', ( $cur_symbol, $cur_id, $ref_seq ) ) . "\n";
				$$info{ $ref_seq }{ 'id' } = $cur_id;
        		$$info{ $ref_seq }{ 'symbol' } = $cur_symbol;       
			}
    	}   
    	elsif( $ref_num lt $cur_num ) { 
			if( $cur_seq ) {
				print STDERR "REV_COMP: " . join( ';', ( $cur_symbol, $cur_id, $cur_seq ) ) . " replaced with " .
                	join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) . "\n";
			} else {      
				print STDERR "DUPS: ". join( ';', ( $cur_symbol, $cur_id, $ref_seq ) ) . 
            		" replaced with " . join( ';', ( $$info{ $ref_seq }{ 'symbol' }, $$info{ $ref_seq }{ 'id' }, $ref_seq ) ) . "\n";
    		}
		}
	}
}
