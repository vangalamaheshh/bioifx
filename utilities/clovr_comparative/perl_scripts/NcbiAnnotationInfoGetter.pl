#!/usr/bin/perl
# Name: NcbiAnnotationInfoGetter.pl
# Desc: Go through all the directories and reads the gen bank files
#	and creates a data struc that has file path and gbk file features
#	writes the data struc to the disk
# 	$$root{$id} = [[filepath, seqLength, seqId], [filepath,seqLength, seqid],...]
# Note: A single taxon id can have multiple refseq files
# Dependencies: Storable
# @author: Mahesh Vangala
##################          LIBRARIES & PRAGMAS            ##################

use strict;
use warnings;
use Storable;

##################            GLOBAL VARIABLES             ##################

my $dir = shift @ARGV || die "Usage: program.pl 'path to top directory that contains ref seq info'\n";
my $root;

##################             SUBS DECLARATION            ##################

sub get_genbank_files ($);
sub print_info ();
sub populate_info ($);

##################              MAIN PROGRAM                #################

get_genbank_files($dir);
store($root, "../binary_files/NcbiAnnotationDataStructure") or die "Error in writing the data structure to disk:\n";
#print_info();
print "Success:\n";
exit(0);

##################              END OF MAIN                 #################

sub get_genbank_files ($) {
	my ($cur_dir) = @_;
	foreach(<$cur_dir/*>) {
		if(-d $_) {
			get_genbank_files($_);
		}
		elsif($_ =~ /\.gbk$/) {
			populate_info($_);
		}
	}
}

sub populate_info ($) {
	open(FH,"<$_") or die "Error in opening the file, $_, $!\n";
  	my ($seq_length, $seq_id);
    while(my $line = <FH>) {
    	if($line =~ /^LOCUS\s+(\S+)\D+(\d+)/) {
    		($seq_length, $seq_id) = ($2, $1);	
    	}
    	elsif($line =~ /^\s+\/db_xref="(.+)"/) {
    		push @{$$root{$1}}, [($_, $seq_length, $seq_id)];
    		close FH;
    		last;
    	}
    }
}

sub print_info () {
	while(my ($key, $value) = each %$root) {
		print $key,"\n";
		foreach(@$value) {
			print join("\t",@{$_}),"\n";	
		}
	}
}

