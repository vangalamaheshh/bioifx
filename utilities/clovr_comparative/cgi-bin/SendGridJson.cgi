#!/usr/bin/perl
# Name: SendGridJson.pl
# Purpose: Based on the node id received from the server, sends the ncbi annotation
#		information back to the client
# Input: Ncbi tree node id
# Output : Ncbi annotation info in JSON format
# Dependencies : Storable, JSON::PP, CGI
# @author : Mahesh Vangala
###################################################################################
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Storable;
use JSON::PP;
use Data::Dumper;

my $CLIENT_ID = 'id';
my $CHILDREN = 'children';
my $TEXT = 'text';
my $ANNOT = 'ANNOTATION_INFO';
my $NAME = 'orgName';
my $LEN = 'seqLen';
my $REF_ID = 'refseqId';
my $SENTINEL = 'userData_';


my $q = new CGI;
my $params = $q->Vars;

my $root;

print "Content-type:text/html\n\n";

if($$params{$CLIENT_ID} =~ /^$SENTINEL(.+)/g) {
	$root = retrieve('/tmp/NcbiUserDataStructure') or
		die "Error retreiving the NcbiUserDataStructure\n";
	my $refInfo = getUserDataAnnotation($1);
	print encode_json({info => $refInfo});
}
else {
	$root = retrieve('/mnt/data/comparative-data-1.0/binary_files/NcbiTreeWidgetDataStructure') or 
		die "Error retrieving the data structure $!\n";
	my $refInfo = getAnnotationInfo($$params{$CLIENT_ID});
	print encode_json({info => $refInfo});
}
exit(0);

sub getUserDataAnnotation {
	my ($id) = @_;
	my $refArray = [];
	if($id eq 'userData') {
		while(my ($key,$value) = each %$root) {
			push @$refArray, @{getUserDataAnnotation($key)};
		}	
	}
	else {
		foreach(@{$$root{$id}{$ANNOT}}) {
			my $refHash = {};
			my @mapName = split(" ", ${$$root{$id}{'MAP_INFO'}}[0]);
			my $mapString = join(" ", @mapName[0..scalar@mapName-2]);
 			$$refHash{$NAME} = $mapString;
			$$refHash{$LEN} = $$_[1];
			$$refHash{$REF_ID} = '<a href=\'http://www.ncbi.nlm.nih.gov/sites/entrez?db=genome&cmd=search&term='.$$_[2].'\' target=\'_blank\'>'.$$_[2].'</a>';
			push @$refArray, $refHash;
		}
	}
	return $refArray;
}

sub getAnnotationInfo {
	my ($id) = @_;
	my $refArray = [];
	while(my ($key,$value) = each %{$$root{$id}{$CHILDREN}}) {
		push @$refArray, @{getAnnotationInfo($key)};
	}
	if($$root{$id}{$ANNOT} && scalar@{$$root{$id}{$ANNOT}} > 0) {
		foreach(@{$$root{$id}{$ANNOT}}) {
			my $refHash = {};
			$$refHash{$NAME} = $$root{$id}{$TEXT};
			$$refHash{$LEN} = $$_[1];
			$$refHash{$REF_ID} = '<a href=\'http://www.ncbi.nlm.nih.gov/sites/entrez?db=genome&cmd=search&term='.$$_[2].'\' target=\'_blank\'>'.$$_[2].'</a>';
			push @$refArray, $refHash;
		}
	}
	return $refArray;
}
