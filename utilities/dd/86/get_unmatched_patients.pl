#!/usr/bin/env perl
use strict;
use warnings;

my $prev_patients_file = shift @ARGV;
my $new_patients_file = shift @ARGV;
my $icd10_file = shift @ARGV;

my $prev = get_info($prev_patients_file);
my $new = get_info($new_patients_file);
my $icd = get_icd($icd10_file);

print_info($prev, $new, $icd);
exit $?;

sub get_info {
  my( $file, $token) = @_;
  open( FH, "<$file" ) or die "Error opening the file, $file, $!\n";
  my $info = {};
  while( my $line = <FH> ) {
    chomp $line;
    if( $line =~ /^(\d+?),\"(.+?)\"/ ) {
      my $mrn = $1;
      my $pat_name = $2;
      $pat_name =~ s/\s+//g;
      $$info{ $mrn } = {'name' => $pat_name, 'line' => $line};
    }
  }
  close FH or die "Error closing the file, $file, $!\n";
  return $info;
}

sub get_icd {
  my( $file ) = @_;
  open( FH, "<$file" ) or die "Error in opening the file, $file $!\n";
  my $info = {};
  my $header = <FH>;
  while( my $line = <FH> ) {
    chomp $line;
    my( $mrn, $icd ) = split(",", $line);
    push @{$$info{$mrn}}, $icd;
  }
  close FH or die "Error in closing the file, $file, $!\n";
  return $info;
}

sub print_info {
  my( $prev, $new, $icd ) = @_;
  foreach my $cur_pat( keys %$new ) {
    if( exists $$prev{$cur_pat} ) {
      print STDERR $cur_pat, ',"', $$new{$cur_pat}{'name'}, '"', ',"', $$prev{$cur_pat}{'name'}, '"', "\n";
    } else {
      print STDOUT $$new{$cur_pat}{'line'}, ',"', join(",", @{$$icd{$cur_pat}}), '"', "\n";
    }
  }
}
