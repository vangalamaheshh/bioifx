#!/usr/bin/env perl

# vim: syntax=perl tabstop=4 expandtab

#---------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Sep, 6, 2016
# --------------------------------

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $options = parse_options();
my $info = get_info( $$options{ 'metasheet' } );
$info = update_file_name( $info );
print_info( $info, $$options{ 'universal_primer' }, $$options{ 'library' } );
exit $?;

sub parse_options {
    my $options = {};
    GetOptions( $options, 'metasheet|m=s', 'universal_primer|u:s', 'library|l:s', 'help|h' );

    my $usage = "$0 
                    Required Parameters:
                    --------------------
                    <--metasheet|-m> - A csv file with 'SampleName, BarCode, IlluminaAdapter'
                    
                    Optional Parameters:
                    --------------------
                    [--universal_primer|-u]
                    [--library|-l]";

    unless( $$options{ 'metasheet' } ) { print STDERR $usage, "\n"; exit 1; }
    if( $$options{ 'help' } ) { print STDERR $usage, "\n"; exit 1; }
    if( not $$options{ 'universal_primer' } ) { $$options{ 'universal_primer' } = 'TCTTGTGGAAAGGACGAAACACCG'; }
    if( not $$options{ 'library' } ) { $$options{ 'library' } = 'Human_A'; }
    return $options;
}

sub get_info {
    my( $file ) = @_;
    my $info = {};
    my $sample_names = {};
    my $ia_bc_combo = {}; # illumina_adapter and barcode combo duplicates check
    open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
    my $header = <FH>; 
    while( my $line = <FH> ) {
        chomp $line;
        my( $sample_name, $barcode, $illumina_adapter ) = split(',', $line );
        $barcode = uc($barcode);
        #$illumina_adapter = uc($illumina_adapter); #rev_comp(uc( $illumina_adapter ));
        $$info{ $illumina_adapter }{ 'barcodes'}{ $barcode }{ 'sample_name' } = $sample_name;
        if( exists $$ia_bc_combo{ $illumina_adapter . $barcode } ) {
            die "ERROR: Duplicate Adapter and Barcode combo: adapter - ${illumina_adapter}; barcode - ${barcode}\n";
        }
        if( exists $$sample_names{ $sample_name } ){
            die "ERROR: Duplicate Sample Name: ${sample_name}\n";
        }
        $$sample_names{ $sample_name }++;
        $$ia_bc_combo{ $illumina_adapter . $barcode }++;
    }
    close FH or die "Error in closing the file, $file, $!\n";
    return $info;
}

sub rev_comp {
    my( $adapter ) = @_;
    $adapter =~ tr/AGCT/TCGA/;
    $adapter = reverse( $adapter );
    return $adapter;
}

sub update_file_name {
    my( $info ) = @_;
    my $data_dir = undef;
    opendir($data_dir, "./data") || die "Error opening directory, data/, $!\n";
    my @files = grep { /^\d+.*.fastq.gz/ && -l "./data/$_" } readdir($data_dir);
    closedir $data_dir;
    my %adapters = map { $_ => 0 } keys %$info;
    foreach my $file( @files ) {
        my( $adapter ) = ( $file =~ /_(\w+?)_EP/ );
        if( exists $adapters{ $adapter } ) {
            $adapters{ $adapter } = 1;
            $$info{ $adapter }{ 'filename' } = $file;           
        }
        else {
            print STDERR "WARN: $file will not be processed since the metasheet doens't contain the adapter, $adapter\n";
        }
    }
    my $NOT_FOUND = 0; #false to begin with
    foreach my $adapter( keys %adapters ) {
        unless( $adapters{ $adapter } ) { 
            print STDERR "ERROR: There is no file with illumina adapter, $adapter\n";
            $NOT_FOUND = 1;
        }
    }
    die "Exiting ....!\n" if( $NOT_FOUND );   
    return $info;
}

sub print_info {
    my( $info, $uni, $lib ) = @_;
    my @header = qw(SampleName FileName BarCode IlluminaAdapter UniversalPrimer Library);
    print STDOUT join(",", @header), "\n";
    foreach my $adapter( keys %$info ) {
        foreach my $barcode( keys %{ $$info{ $adapter }{ 'barcodes' } } ) {
            print STDOUT join(',', ( $$info{ $adapter }{ 'barcodes' }{ $barcode }{ 'sample_name' },
                                    $$info{ $adapter }{ 'filename' }, $barcode, $adapter, $uni, $lib )), "\n";           
        }
    }   
}
