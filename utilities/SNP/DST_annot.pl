#!/usr/bin/env perl
use strict;
use warnings;

my $csv_info_file = $ARGV[0];
my $variant_db_file = $ARGV[1];

my $variant_db_info = get_db_info($variant_db_file);
process_data($csv_info_file, $variant_db_info);
exit $?;

sub get_db_info {
  my($file) = @_;
  my $info = {};
  open(FH, "<$file") or die "Error opening the file, $file, $!\n";
  my $header = <FH>;
  while(my $line = <FH>) {
    chomp $line;
    my @array = split(",", $line);
    $$info{$array[0]} = $line;
  }
  close FH or die "Error in clsoing the file, $file, $!\n";
  return $info;
}

sub process_data {
  my($file, $db) = @_;
  open(FH, "<$file") or die "Error in openign the file, $file, $!\n";
  my $header = <FH>;
  while(my $line = <FH>) {
    chomp $line;
    my($num) = ($line =~ /^(\d+)/);
    my($cosm) = ($line =~ /(COSM\d+)/);
    my($rs) = ($line =~ /(rs\d+)/);
    my $info = [$num];
    if($cosm or $rs) {
      if($cosm and exists $$db{$cosm}) {
        push @$info, $$db{$cosm};
      }
      if($rs and exists $$db{$rs}) {
        push @$info, $$db{$rs};
      }
    }
    print STDOUT join(",", @$info), "\n";
  }
  close FH or die "Error in closing the file, $file, $!\n";
}
