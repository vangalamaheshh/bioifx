#!/usr/bin/perl
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use JSON::PP;
use Storable;
use Data::Dumper;

my $QUERY = 'query';
my $NODE = 'selectedNode';

my $root = retrieve('/mnt/data/comparative-data-1.0/binary_files/AutoSearchDataStructure')
	|| die "Error in retrieving the data structure, $!\n";
	
my $q = new CGI;
my $params = $q->Vars;

print "Content-type: text/html\n\n";
if($$params{$NODE}) {
	print encode_json($$root{$$params{$NODE}});
}
elsif($$params{$QUERY}) {
	my $nodeNamesList = getMatchedNamesList($root);
	print encode_json($nodeNamesList);	
}
else {
	my $nodeNamesList = getNodeNamesList($root);
	print encode_json($nodeNamesList);
}
exit(0);

sub getNodeNamesList {
	my ($info) = @_;
	my $ref;
	foreach(sort keys %$info) {
		push @$ref, [$_];	
	}
	return $ref;
}

sub getMatchedNamesList {
	my ($info) = @_;
	my $ref;
	foreach(sort keys %$info) {
		if($_ =~ /$$params{$QUERY}/gi) {
			push @$ref, [$_];
		}
	}
	return $ref;
}


