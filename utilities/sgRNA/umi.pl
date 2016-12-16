#!/usr/bin/env perl
# vim: syntax=perl tabstop=4 expandtab

use strict;
use warnings;
use Getopt::Long;

#--------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 26, 2016
#--------------------

my $SPACER_1 = 'GTGATTGCTTGTGACGCCTT';
my $SPACER_2 = 'CAGGAGTGAGATGACAGGAG';
my $sgRNA_FILE = '/zfs/cores/mbcf/mbcf-storage/devel/umv/sandor_sgRNA_1018.csv';

my $options = parse_options();
my $sgInfo = get_sgRNA_info( $sgRNA_FILE );
my $info = get_info( $$options{ 'fastq_file' }, $sgInfo );
print_info( $info, $sgInfo );

sub parse_options {
	my $options = {};
	GetOptions( $options, '--fastq_file|f=s@', 'help|h' );
	unless( $$options{ 'fastq_file' } ) {
        my $usage = "$0 <--fastq_file|-f> [--fastq_file|-f, ...]";
        print STDERR $usage,"\n";		
	}
    return $options;
}

sub get_sgRNA_info {
    my( $sg_file ) = @_;
    my $sg_info = {};
    open( FH, "<$sg_file" ) or die "Error in opening the file, $sg_file, $!\n";
    my $header = <FH>;
    while( my $line = <FH> ) {
        chomp $line;
        my( $name, $seq ) = split( ",", $line );
        $$sg_info{ $seq } = $name;
    }
    close FH or die "Error in closing the file, $sg_file, $!\n";
    return $sg_info;
}

sub get_info {
    my( $fastq_list, $sgInfo ) = @_;
    my $info = {};
    foreach my $fastq_file( @$fastq_list ){
        my( $sample_name ) = ( $fastq_file =~ /(.+)_R\d_/ );
        print STDERR "INFO: Processing $sample_name\n";
        open( FH, "zcat $fastq_file |" ) or die "Error opening gzipped file, $fastq_file, $!\n";
        while( my $header = <FH> ) {
            my $seq = <FH>;
            my $sep = <FH>;
            my $ascii = <FH>;
            if( $seq =~ /(........)${SPACER_1}(........)${SPACER_2}/ ) {
                my $barcode = $1.$2;
                $$info{ 'barcode_combo' }{ $barcode }{ 'count' }++;
                foreach my $sgRNA( keys %$sgInfo ) {
                    my $rev_comp = $sgRNA;
                    $rev_comp =~ tr/AGCT/TCGA/;
                    $rev_comp = reverse($rev_comp);
                    if( $seq =~ /($sgRNA)/ or $seq =~ /($rev_comp)/ ) {
                        $$info{ 'barcode_combo' }{ $barcode }{ 'sgRNA' }{ $1 }++;
                    }
                }
            }
        }
        close FH or die "Error closing file, $fastq_file, $!\n";
    }
    return $info;
}

sub print_info {
    my( $info, $sgInfo ) = @_;
    my @header = qw(BarcodeCombo Count Unique_sgRNA sgRNA);
    print STDOUT join( ",", @header ), "\n"; 
    foreach my $barcode( keys %{ $$info{ 'barcode_combo' } } ) {
        my @sgRNAs = ();
        my @sgCounts = ();
        foreach my $sgRNA( keys %{ $$info{'barcode_combo'}{$barcode}{'sgRNA'} } ){
            push @sgRNAs, $$sgInfo{ $sgRNA };
            push @sgCounts, $sgRNA . '(' . $$info{ 'barcode_combo' }{ $barcode }{ 'sgRNA' }{ $sgRNA } . ')';
        }
        my $sgRNA_str = join(';', @sgRNAs);
        my $sgcount_str = join(';', @sgCounts);
        print STDOUT join( ",", ( $barcode, $$info{ 'barcode_combo' }{ $barcode }{ 'count' },
            $sgRNA_str, $sgcount_str ) ), "\n";
    }
}
