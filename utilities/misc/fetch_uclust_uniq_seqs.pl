#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $options = parse_options( );
my $file_info = get_file_info( $$options{ 'uclust_dir' } );

my $fasta_token = '.filtered.sorted.fasta';
my $uclust_token = '.uclust.out';
my $out_token = '.uclust.uniq.fasta';

foreach my $file( keys %$file_info ) {
    print STDERR "INFO: Processing $file\n";
    my $uclust_uniq_headers = get_uclust_uniq_headers( $$options{ 'uclust_dir' } . '/' . $file . $uclust_token );
    print STDERR "INFO: Number of unique sequences found: ", scalar keys %$uclust_uniq_headers, "\n";
    my $in_file = $$options{ 'uclust_dir' } . '/' . $file . $fasta_token;
    my $out_file = $$options{ 'uclust_dir' } . '/' . $file . $out_token;
    open( IN, "<$in_file" ) or die "Error in opening the file, $in_file, $!\n";
    open( OUT, ">$out_file" ) or die "Error in writing to the file, $out_file, $!\n";
    my $counter = 0;
    while( my $header = <IN> ) {
        unless( $header =~ /^>/ ) {
            die "Error: $header doesn't start with token '>'\n";
        }
        my $seq = <IN>;
        chomp $header;
        chomp $seq;
        if( exists $$uclust_uniq_headers{ $header } ) {
            print OUT $header, "\n";
            print OUT $seq, "\n";
            $counter++;
        }
    }
    close IN or die "Error in closing the file, $in_file, $!\n";
    close OUT or die "Error in closing the file, $out_file, $!\n";
    unless( scalar keys %$uclust_uniq_headers == $counter ) {
        print STDERR "Error: Counter: $counter is not equal to UCLUST_UNIQ_HEADERS: ", scalar keys %$uclust_uniq_headers, "\n";
    } 
    print STDERR "INFO: Done\n";
}

exit $?;

sub parse_options {
    my $options = {};
    GetOptions( $options, 'uclust_dir|u=s', 'help|h' );
    unless( $$options{ 'uclust_dir' } ) {
        die "Usage: $0 <--uclust_dir|-u>\n";
    }
    return $options;
}

sub get_file_info {
    my( $dir ) = @_;
    my $info = {};
    opendir( DIR, $dir ) or die "Error in opening directory, $dir, $!\n";
    while( my $file = readdir( DIR ) ) {
        chomp $file;
        if( $file =~ /(.+?)\./ ) {
            next if( $file =~ /^\./ );
            if( $1 ) {
                $$info{ $1 }++;
            }
        }
    }
    closedir DIR or die "Error in closing directory, $dir, $!\n";
    return $info;
}

sub print_file_info {
    my( $file_info ) = @_;
    foreach my $file_name( keys %$file_info ) {
        print $file_name, "\n";
    }
}

sub get_uclust_uniq_headers {
    my( $uclust_file ) = @_;
    my $uclust_uniq_headers = {};
    open( UCLUST, "<$uclust_file" ) or die "Error in opening the file, $uclust_file, $!\n";
    while( my $line = <UCLUST> ) {
        next if $line =~ /^#/;
        chomp $line;
        my @array = split( "\t", $line );
        if( $array[ 0 ] eq 'S' ) {
            $$uclust_uniq_headers{ '>' . $array[ 8 ] }++;
        }
    }
    close UCLUST or die "Error in closing the file, $uclust_file, $!\n";
    return $uclust_uniq_headers;
}
