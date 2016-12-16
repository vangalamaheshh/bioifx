#!/usr/bin/perl
# Name: UserDataAnnotationInfoGetter.cgi
# IMPORTANT: This program logic is similar to NcbiAnnotationInfoGetter.pl
#		The only reason, why I created this script is, we have to read the
#		user data, get the refseq info and map info all at once and write 
#		to the disc. NcbiAnnotationInfoGetter.pl only adds ref seq info and
#		the map info is added in NcbiTaxaAnnotation.pl. 
# Desc: Go through all the directories and reads the gen bank files
#	and creates a data struc that has file path and gbk file features
#	writes the data struc to the disk
# 	$$root{$id}{$ANNOT} = [[filepath, seqLength, seqId], [filepath,seqLength, seqid],...]
#	$$root{$id}{$MAP} = [org1, org2, ...]
# Note: A single taxon id can have multiple refseq files
# Dependencies: Storable
# @author: Mahesh Vangala
##################          LIBRARIES & PRAGMAS            ##################

use strict;
use warnings;
use Storable;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use JSON::PP;

##################               CONSTANTS                 ##################

my $ANNOT = 'ANNOTATION_INFO';
my $MAP = 'MAP_INFO';

##################            GLOBAL VARIABLES             ##################

#my $dir = shift @ARGV || die "Usage: program.pl 'path to top directory that contains ref seq info'\n";
my $root;

##################             SUBS DECLARATION            ##################

sub get_genbank_files ($);
sub print_info ();
sub populate_info ($);

##################              MAIN PROGRAM                #################

my $q = new CGI;
my $params = $q->Vars;
print "Content-type: text/html\n\n";

my $fileListRef = decode_json($$params{'fileList'});
foreach(@$fileListRef) {
	populate_info($_);
}
#get_genbank_files($dir);
store($root, "/tmp/NcbiUserDataStructure") or die "Error in writing the data structure to disk:\n";
#print_info();
print encode_json([JSON::PP::true]);
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
	return unless(-e $_);
	open(FH,"<$_") or die "Error in opening the file, $_, $!\n";
  	my ($seq_length, $seq_id, $mapOrgName);
    while(my $line = <FH>) {
    	if($line =~ /^LOCUS\s+(\S+)\D+(\d+)/) {
    		($seq_length, $seq_id) = ($2, $1);	
    	}
	elsif($line =~ /^\s+ORGANISM\s+(.+)/) {
	     	my @array = split(" ", $1);
                $mapOrgName = join(" ", (@array,$array[scalar@array-1]));
	}
    	elsif($line =~ /^\s+\/db_xref="(.+)"/) {
    		push @{$$root{$1}{$ANNOT}}, [($_, $seq_length, $seq_id)];
			push @{$$root{$1}{$MAP}}, $mapOrgName;
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

