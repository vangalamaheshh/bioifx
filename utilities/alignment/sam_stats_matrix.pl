#!/usr/bin/env perl
# vim: syntax=perl tabstop=4 expandtab

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my $sample_list = [map { my($temp) = (basename($_) =~ /(.+?)\./); $temp; } @{$$options{'sample'}}];
my $info = get_matrix($$options{'sample'});
print_info($info, $sample_list);
exit(0);

sub parse_options {
    my $options = {};
    GetOptions($options, 'sample|s=s@', 'help|h');
    unless($$options{'sample'}) {
        print STDERR "Usage: $0 <--sample|-s> [--sample|-s]\n";
        exit(1);
    } 
    return $options;
}

sub get_matrix {
    my($ref_samples) = @_;
    my $info = {};
    foreach my $sample_file(@$ref_samples) {
        my ($sample) = (basename($sample_file) =~ /(.+?)\./);
        open(FH, "<$sample_file") or die "Error opening file, $sample_file, $!\n";
        my($total, $unique) = (undef, undef);
        while(my $line = <FH>) {
            chomp $line;
            if($line =~ /raw\stotal\ssequences:\s(\d+)/) {
                $total = $1;
            }
            elsif($line =~ /reads\smapped:\s(\d+)/) {
                $unique = $1;
            }
            last if $total and $unique;
        }
        close FH or die "Error closing file, $sample_file, $!\n";
        $$info{$sample}{'Total'} = $total;
        $$info{$sample}{'Unique'} = $unique;
    }
    return $info;
}

sub print_info {
    my($matrix, $samples) = @_;
    my $header = ",TotalReadCount,UniqueReadCount\n";
    print STDOUT $header;
    foreach my $sample(@$samples) {
        print STDOUT join(",", ($sample, $$matrix{$sample}{'Total'},
                $$matrix{$sample}{'Unique'})), "\n";
    }
}
