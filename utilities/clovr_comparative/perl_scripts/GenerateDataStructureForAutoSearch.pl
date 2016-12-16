#!/usr/bin/perl
use strict;
use warnings;
use Storable;

my $TEXT = 'text';
my $ROOT = "/";
my $CHILDREN = 'children';
my $LEAF = 'leaf';
my $name = '';
my $finalRef;

my $root = retrieve("../binary_files/NcbiTreeWidgetDataStructure") or
	die "Error in retrieving the data structure, $!\n";
	
getFinalDataStructure($root);
#printInfo($finalRef);
store($finalRef, '../binary_files/AutoSearchDataStructure') or
	die "Error in writing the data structure to the disc, $!\n";
print "Success\n";
exit(0);

###############

sub getFinalDataStructure {
	my ($root) = @_;
	helperSub($$root{$ROOT}{$CHILDREN}, ());
}

sub helperSub {
	my ($info, @lineage) = @_;
	while(my ($key,$value) = each %$info) {
		unless($$root{$key}{$LEAF}) {
			helperSub($$root{$key}{$CHILDREN}, (@lineage,$key)); 
		}
		##################     IMPORTANT        ###################
		# Now we should append the current key as the last node id
		# so that the node typed by the user will also expand
		###########################################################
		$$finalRef{$$root{$key}{$TEXT}} = [@lineage, $key];		
	}
}

sub printInfo {
	my ($info) = @_;
	while(my ($key,$value) = each %$info) {
		print $key,"\t",join(", ",@$value),"\n";
	}
}
