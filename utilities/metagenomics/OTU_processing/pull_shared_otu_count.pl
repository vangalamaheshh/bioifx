#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

=head1 SYNOPSIS

	<programname.pl> --raw_otu_file|-r --sample_names_file|-s 
			 --otu_taxa_file|-t --output_dir|-o [--help|-h]

	PARAMETERS:
		
		REQUIRED:
		--raw_otu_file|-r
			Raw OTU file from uclust_picked_otus component in Qiime.

		--sample_names_file|-s
			File with sample names

		--otu_taxa_file|-t
			Output file from RDP_Classifier

		--output_dir|-o
			Output directory name. The output files are written to this dir.

		OPTIONAL:
		--help|-h
			Prints this message and exit.

=head1 AUTHOR
	Mahesh Vangala
	vangalamaheshh@gmail.com
=cut

#---------------------------------------------------------------------#
# MAIN PROGRAM
#---------------------------------------------------------------------#
my $options = parse_options( );
my $otu_taxa_info = get_otu_taxa_info( $$options{ 'otu_taxa_file' } );
my $info = initialize_hash( $$options{ 'sample_names_file' } );
$info = process( $info, $$options{ 'raw_otu_file' }, $otu_taxa_info );
print_to_files( $$options{ 'output_dir' } );
exit $?;
#---------------------------------------------------------------------#


sub print_to_files {
	my( $output_dir ) = @_;
	unless( -d $output_dir ) {
		mkdir $output_dir or die "Cannot create directory, $output_dir, $!\n";
	}
	my $stats_file = "$output_dir/stats.txt";
	print_stats( $info, $stats_file );
	print_taxa( $$info{ 'root' }, $output_dir, 'root' );
}

sub print_taxa {
	my( $node, $output_dir, $node_name ) = @_;
	if( defined $$node{ 'children' } ) {
		my $output_file = "$output_dir/$node_name.shared.taxa";
		open( OFH, ">$output_file" ) or die "Error in writing to the file, $output_file, $!\n";
		foreach my $index( 0 .. scalar @{ $$node{ 'otus' } } - 1 ) {
			print OFH join( "\t", ( ${ $$node{ 'otus' } }[ $index ], ${ $$node{ 'taxa_info' } }[ $index ] ) ), "\n"; 
		}
		close OFH or die "Error in closing the file, $output_file, $!\n";
		foreach my $child( keys %{ $$node{ 'children' } } ) {
			print_taxa( $$node{ 'children' }{ $child }, $output_dir, $child );
		}
	}
}

sub print_stats {
	my( $info, $stats_file ) = @_;
	open( OFH, ">$stats_file" ) or die "Error in writing to the file, $stats_file, $!\n";
	print OFH "root", "\t", scalar @{ $$info{ 'root' }{ 'otus' } }, "\n";
	print_info( $$info{ 'root' }{ 'children' } );
	close OFH or die "Error in closing the file, $stats_file, $!\n";
}

sub get_otu_taxa_info {
	my( $otu_taxa_file ) = @_;
	my $info = {};
	open( FH, "<$otu_taxa_file" ) or die "Error in opening the file, $otu_taxa_file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my( $otu, $read_id, $taxa, $val ) = ( $line =~ /(.+?)\s+(.+?)\s+(.+?)\s+(.+)/ );
		$$info{ $otu } = $taxa;
	}
	close FH or die "Error in clsing the file, $otu_taxa_file, $!\n";
	return $info;
}

sub print_info {
	my( $root ) = @_;
	while( my( $name, $ref_hash ) = each %$root ) {
		if( defined $$ref_hash{ 'count' } ) {
			print OFH $name, "\t", scalar @{ $$ref_hash{ 'otus' } }, "\n";
			print_info( $$ref_hash{'children'} );
		}
	}
}

sub process {
	my( $info, $otu_file, $otu_taxa_info ) = @_;
	open( FH, "<$otu_file" ) or die "Error opening the file, $otu_file, $!\n";
	while( my $line = <FH> ) {
		chomp $line;
		my( $otu_id, @otu_list ) = split( "\t", $line );
		get_shared_otu_count( $$info{ 'root' }, undef, $otu_id, \@otu_list, $otu_taxa_info );
		initialize_count_to_zero( $$info{ 'root' } );
	}
	close FH or die "Error clsong the file, $otu_file, $!\n";
	return $info;
}

sub initialize_count_to_zero {
	my( $root ) = @_;
	if( defined $$root{ 'count' } ) {
		$$root{ 'count' } = 0;
		foreach my $child( keys %{ $$root{ 'children' } } ) {
			initialize_count_to_zero( $$root{ 'children' }{ $child } );
		}
	}
}

sub get_shared_otu_count {
	my( $parent, $key_name, $otu_id, $otu_list, $otu_taxa_info ) = @_;
	if( defined $$parent{ 'children' } ) {
		foreach my $child( keys %{ $$parent{ 'children' } } ) {
			$$parent{ 'count' } += get_shared_otu_count( $$parent{ 'children' }{ $child }, $child, $otu_id, $otu_list, $otu_taxa_info );
		}
	}
	else {
		if( join( "\t", @$otu_list ) =~ /$key_name/ ) {
			return 1;
		}
		else {
			return 0;
		}
	}

	if( $$parent{ 'count' } == scalar keys %{ $$parent{ 'children' } } ) {
		push @{ $$parent{ 'otus' } }, $otu_id;
		push @{ $$parent{ 'taxa_info' } }, $$otu_taxa_info{ $otu_id };
		return 1;
	}
	return 0;
}

sub initialize_hash {
	my( $file ) = @_;
	open( FH, "<$file" ) or die "Error opening the file, $file, $!\n";
	my $info = { root => {
				count => 0,
				otus => [],
				taxa_info => [],
				children => {}
				} 
			};
	while( my $line = <FH> ) {
		chomp $line;
		my( $lake, $time_point ) = ( $line =~ /(\w+)_\w(\d)/ );
		unless( defined $$info{ 'root' }{ 'children' }{ $lake } ) {
			$$info{ 'root' }{ 'children' }{ $lake } = { count => 0, otus => [], taxa_info => [], children => {} };
		}
		unless( defined $$info{ 'root' }{ 'children' }{ $lake }{ 'children' }{ $lake.$time_point } ) {
			$$info{ 'root' }{ 'children' }{ $lake }{ 'children' }{ $lake.$time_point } = { count => 0, otus => [], taxa_info => [], children => {} };
		}
		$$info{ 'root' }{ 'children' }{ $lake }{ 'children' }{ $lake.$time_point }{ 'children' }{ $line } = {};
	}
	close FH or die "Error clsing the file, $file, $!\n";
	return $info;
}

sub parse_options {
	my $options = {};
	GetOptions( $options, 'raw_otu_file|r=s', 'sample_names_file|s=s', 'otu_taxa_file|t=s', 'output_dir|o=s', 'help|h' ) 
					|| pod2usage( { -output => \*STDERR, -exitval => 1, -verbose => 1 } );
	if( $$options{ 'help' } or not $$options{ 'sample_names_file' } or not $$options{ 'raw_otu_file' } 
		or not $$options{ 'otu_taxa_file' } or not $$options{ 'output_dir' } ) {
		pod2usage( { -output => \*STDERR, -exitval => 1, -verbose => 1 } );
	}
	return $options;
}

