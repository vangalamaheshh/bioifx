#!/usr/bin/env perl

# vim: syntax=perl tabstop=4 expandtab

#----------------------
#@author: Mahesh Vangala
#@email: vangalamaheshh@gmail.com
#@date: Aug, 05, 2016
#----------------------

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my ($matrix) = get_matrix( $$options{ 'fastqc_file' } );
print_matrix($matrix);
exit $?;

sub parse_options {
    my $options = {};
    GetOptions( $options, 'fastqc_file|f=s@', 'help|h' );
    unless( $$options{ 'fastqc_file' } ) {
        my $usage = "$0 <--fastqc_file|-f> [--fastqc-file|-f ...]";
        print STDERR $usage, "\n";
        exit 1;
    }	
    return $options;
}

sub get_matrix {
    my( $file_list ) = @_;
    my $info = {};
    foreach my $file( @$file_list ){
        my $sample = basename($file);
        $sample = substr($sample, 0, -4);
        open( FH, "unzip -c $file ${sample}/summary.txt |" ) or
            die "Error executing command, $!\n";
        my $line = <FH>;
        $line = <FH>;
        while( $line = <FH> ) {
            chomp $line;
            next if $line =~ /^\s*$/;
            my( $flag, $type ) = split( "\t", $line );
            $type =~ s/\s+/_/g;
            $$info{ $type }{ $sample } = $flag;
        }     
        close FH or die "Error closing the pipe, $!\n";
        open( FH, "unzip -c $file ${sample}/fastqc_data.txt |" ) or
            die "Error executing command, $!\n";
        while( $line = <FH> ) {
            next if $line !~ /Total/;
            chomp $line;
            $line =~ s/#//g;
            my( $type, $num ) = split( "\t", $line );
            $type =~ s/\s+/_/g;
            $$info{ $type }{ $sample } = $num;
        }
        close FH or die "Error closing the pipe, $!\n";
    }
    return $info; 
}

sub print_matrix {
    my( $info ) = @_;
    my @samples = keys %{$$info{'Basic_Statistics'}};
    my @short_names = @samples;
    foreach my $name(@short_names) { $name =~ s/((\d|_)+)_[^\d]//; $name =~ s/_\d+_fastqc//; }
    print STDOUT join(',', ("Feature",@short_names)) . "\n";
    foreach my $feature( keys %$info ) {
        my @str = ();
        foreach my $sample( @samples ) {
            push @str, $$info{ $feature }{ $sample };
        }
        print STDOUT join(',', ( $feature, @str )), "\n";
    }   
}
