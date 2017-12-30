#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift @ARGV;

my $report = get_report($file);
print_report($report);
exit $?;

sub get_report {
  my($file) = @_;
  my $report = {};
  open(FH, "<$file") or die "Error opening the file, $file, $!\n";
  my $line = <FH>;
  OUTERLOOP:
  while($line) {
    if(length((split(",", $line))[0]) == 3) {
      #start new day
      $$report{'days'} += 1;
      INNERLOOP:
      while($line = <FH>) {
        if(length((split(",", $line))[0]) == 3) {
          redo OUTERLOOP;
        }
        $line =~ s/^\s+//;
        if($line =~ /^[1-9]/) {
          #start new activity
          my($topic, $effort) = ($line =~ /^\d+.\s*(.+):\s*(\d+)/);
          if($topic =~ /workshop|seminar|conference|webinar|training/i) {
            $topic = "Workshops";
          }
          elsif($topic =~ /doyle|mamba|philips|mghpcc/i) {
            $topic = "MAMBA";
          }
          elsif($topic =~ /holiday/i) {
            $topic = "Holidays";
          }
          elsif($topic =~ /meeting/i) {
            $topic = "Meetings";
          }
          elsif($topic =~ /presentation/i) {
            $topic = "Presentations";
          }
          elsif($topic =~ /i2b2/i) {
            $topic = "i2b2-tranSMART-shrine";
          }
          elsif($topic =~ /trinet/i) {
            $topic = "TriNetX";
          }
          elsif($topic =~ /gcp|synergist|variant\s+calling/i) {
            $topic = "Synergist";
          }
          $$report{'topic'}{$topic}{'effort'} += $effort;
          $$report{'topic'}{$topic}{'days'} += 1;
        }
      }
    }
  }
  close FH or die "Error closing the file, $file, $!\n";
  return $report;
}

sub print_report {
  my($report) = @_;
  my @header = qw(Topic Effort DaysWorked YearlyPercentage);
  print join(",", @header), "\n";
  foreach my $topic(keys %{$$report{'topic'}}) {
    my @vals = ($topic,
                $$report{'topic'}{$topic}{'effort'},
                $$report{'topic'}{$topic}{'days'});
    push @vals, $vals[1]/$$report{'days'};
    print join(",", @vals), "\n";
  }  
}
