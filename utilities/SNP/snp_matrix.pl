#!/usr/bin/env perl
#vim: syntax=perl tabstop=2 expandtab

#---------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Jan, 11, 2017
#---------------------

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $options = parse_options();
my ($vcf_info, $chrom_list) = get_vcf_info($$options{'variant_file'});
print_report($vcf_info, $$options{'report'});
print_summary_report($vcf_info, $chrom_list, $$options{'summary_report'});
exit $?;

sub parse_options {
  my $options = {};
  GetOptions($options, 'variant_file|v=s@', 'report|r=s', 'summary_report|s=s', 'help|h');
  unless( $$options{'variant_file'} and $$options{'report'} and $$options{'summary_report'}) {
    print STDERR "Usage: $0 <--variant_file|-v> <--report|-r> <--summary_report|-s>\n
                                  [multiple variant files can be given as -v <file1> -v <file2>\n";
    exit(1);
  }
  return $options;
}

sub get_vcf_info {
  my ($vcf_file_list) = @_;
  my $vcf_info = {};
  my $chroms = {};
  foreach my $vcf_file(@$vcf_file_list) {
    my $token = basename($vcf_file);
    $token =~ s/\..+//g;
    open(FH, "<$vcf_file") or die "Error in opening the file, $vcf_file, $!\n";
    while(my $line = <FH>) {
      chomp $line;
      next if($line =~ /^#/);
      my($chrom, $pos, $id, $ref, $alt, $qual, $filter) = split("\t", $line);
      if ($filter eq 'PASS') {
        $$vcf_info{$token}{$chrom}{$pos} = {'id' => $id, 'ref' => $ref, 'alt' => $alt};
        $$chroms{$chrom}++;
      }
    }
    close FH or die "Error closing the file, $vcf_file, $!\n";
  }
  my @chrom_list = sort {$a cmp $b} keys %$chroms;
  return $vcf_info, \@chrom_list;
}

sub print_report {
  my ($info, $file) = @_;
  open(OFH, ">$file") or die "Error in opening the file, $file, $!\n";
  my @header = qw(SampleName Chromosome Position Reference_Allele Alternate_Allele dbSNP_RS);
  print OFH join(",", @header), "\n";
  foreach my $sample(sort {$a cmp $b} keys %$info) {
    foreach my $chrom(sort {$a cmp $b} keys %{$$info{$sample}}) {
      foreach my $pos(sort {$a <=> $b} keys %{$$info{$sample}{$chrom}}) {
        print OFH join(",", ($sample, $chrom, $pos,
          $$info{$sample}{$chrom}{$pos}{'ref'}, $$info{$sample}{$chrom}{$pos}{'alt'},
          $$info{$sample}{$chrom}{$pos}{'id'})), "\n";
      }
    }
  }
  close OFH or die "Error closing the file, $file, $!\n";
}

sub print_summary_report {
  my($info, $chrom_list, $file) = @_;
  open(OFH, ">$file") or die "Error in opening the file, $file, $!\n";
  print OFH join(",", ("SampleName", @$chrom_list)), "\n";
  foreach my $sample (sort {$a cmp $b} keys %$info) {
    my @vals = ();
    foreach my $chrom (@$chrom_list) {
      push @vals, scalar keys %{$$info{$sample}{$chrom}} || 0;
    }
    print OFH join(",", ($sample, @vals)), "\n";
  } 
  close OFH or die "Error in closing the file, $file, $!\n";
}


